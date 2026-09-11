#!/usr/bin/env python3
"""Project saved global context-similarity maps using existing FSL transforms.

Run with a Python containing numpy/nibabel (e.g. $FSLDIR/bin/python):
  python scripts/project_global_searchlight_to_T1.py --subjects 2 3 4 5 6
No registration estimation or statistical inference is performed. Existing
outputs are never overwritten. Nearest-neighbour sampling preserves discrete
support, p values, and native functional-grid effect values.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

import nibabel as nib
import numpy as np

CONTEXTS = ('person', 'food', 'location')
SUFFIXES = ('delta_z', 'null_z', 'p_uncorrected', 'p_fwe_maxstat', 'sig_fwe_p05')
SUPPORT = ('valid_center_mask', 'final_restriction_mask', 'roi_union', 'neighborhood_size')


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def project_subject(subject, root, flirt):
    base = root / 'MRI' / f'subj_{subject}' / 'nifti'
    source = base / 'olf_context_similarity_searchlight' / 'global'
    destination = source / 'native_T1'
    reference = base / 'anat' / 'betT1brain.nii.gz'
    transform = base / 'coreg' / 'func2T1mat'
    prefix = f'olf_context_similarity_subj{subject}'
    keys = list(SUPPORT) + [f'{c}_{s}' for c in CONTEXTS for s in SUFFIXES]
    inputs = {key: source / f'{prefix}_{key}.nii' for key in keys}
    for path in [reference, transform, *inputs.values()]:
        if not path.is_file():
            raise FileNotFoundError(path)
    if destination.exists():
        raise FileExistsError(f'Preserving existing output: {destination}')
    ref = nib.load(reference)
    grid = nib.load(inputs['valid_center_mask'])
    for path in inputs.values():
        image = nib.load(path)
        assert image.shape == grid.shape and np.allclose(image.affine, grid.affine, atol=1e-4), path
    matrix = np.loadtxt(transform)
    assert matrix.shape == (4, 4) and np.isfinite(matrix).all()
    assert abs(np.linalg.det(matrix[:3, :3])) > 1e-8
    commands = []
    records = []
    env = dict(os.environ, FSLOUTPUTTYPE='NIFTI', OMP_NUM_THREADS='1')
    # Publish the folder only after all images and checks have completed.
    with tempfile.TemporaryDirectory(prefix=f'.native_T1_subj{subject}_', dir=source) as workstr:
        work = Path(workstr)
        publish = work / 'output'
        publish.mkdir()
        def save(values, path):
            header = ref.header.copy()
            header.set_data_dtype(np.float32)
            nib.save(nib.Nifti1Image(values.astype(np.float32), ref.affine, header), path)
        valid = None
        for key, path in inputs.items():
            native = nib.load(path)
            values = native.get_fdata(dtype=np.float32)
            is_p = '_p_uncorrected' in key or '_p_fwe_maxstat' in key
            fill = 1 if is_p else 0
            cleaned = np.where(np.isfinite(values), values, fill)
            clean_path = work / 'input.nii'
            clean_header = native.header.copy()
            clean_header.set_data_dtype(np.float32)
            nib.save(nib.Nifti1Image(cleaned, native.affine, clean_header), clean_path)
            sampled_path = work / 'sampled.nii'
            command = [str(flirt), '-in', str(clean_path), '-ref', str(reference),
                       '-applyxfm', '-init', str(transform), '-interp', 'nearestneighbour',
                       '-out', str(sampled_path)]
            subprocess.run(command, check=True, env=env, capture_output=True, text=True)
            commands.append(dict(source=str(path), interpolation='nearestneighbour'))
            image = nib.load(sampled_path)
            assert image.shape == ref.shape and np.allclose(image.affine, ref.affine, atol=1e-4)
            projected = image.get_fdata(dtype=np.float32)
            assert np.isfinite(projected).all()
            if key == 'valid_center_mask':
                valid = projected > 0
                assert valid.any()
            if key not in SUPPORT[:3]:
                projected[~valid] = fill
            if key.endswith('mask') or key == 'roi_union' or key.endswith('sig_fwe_p05'):
                assert np.isin(projected, [0, 1]).all()
            if is_p:
                assert np.all((projected >= 0) & (projected <= 1))
            # With nearest neighbour every target value must be a source value
            # or the declared background value (no invented statistic values).
            assert np.isin(np.unique(projected), np.append(np.unique(cleaned), fill)).all()
            target = publish / f'{prefix}_{key}_T1.nii.gz'
            save(projected, target)
            records.append(dict(source=str(path), source_sha256=sha256(path), output=target.name,
                                background=fill, source_nonzero=int(np.count_nonzero(cleaned)),
                                target_nonzero=int(np.count_nonzero(projected))))
        for context in CONTEXTS:
            def read(suffix):
                return nib.load(publish / f'{prefix}_{context}_{suffix}_T1.nii.gz').get_fdata(dtype=np.float32)
            sig = read('sig_fwe_p05') > 0
            p = read('p_fwe_maxstat')
            assert np.array_equal(sig, valid & (p < .05))
            native_sig = nib.load(inputs[f'{context}_sig_fwe_p05']).get_fdata() > 0
            assert not native_sig.any() or sig.any(), 'Significance support disappeared'
            for statistic in ('delta_z', 'null_z'):
                save(np.where(sig, read(statistic), 0), publish / f'{prefix}_{context}_{statistic}_sig_fwe_p05_T1.nii.gz')
        metadata = dict(subject=subject, created_utc=datetime.now(timezone.utc).isoformat(),
                        producing_script=str(Path(__file__).resolve()), script_sha256=sha256(Path(__file__)),
                        reference=str(reference), reference_sha256=sha256(reference),
                        transform=str(transform), transform_sha256=sha256(transform),
                        transform_matrix=matrix.tolist(), interpolation='nearestneighbour',
                        valid_T1_centers=int(valid.sum()), maps=records,
                        provenance='Source MAT results remain in parent directory; inference not recomputed.',
                        checks='T1 geometry, finite values, binary masks, source-value preservation, p/FWE-mask agreement passed.')
        (publish / 'projection_manifest.json').write_text(json.dumps(metadata, indent=2) + '\n')
        (publish / 'README.md').write_text(
            f'# Subject {subject}: global searchlight in native T1 space\n\n'
            f'Reference: `{reference}`\n\nTransform: `{transform}`\n\n'
            'All saved functional maps were resampled with FSL FLIRT using the existing '
            'functional-to-T1 affine and nearest-neighbour interpolation. No registration '
            'was estimated and no inference was rerun. Source maps are preserved.\n\n'
            'Effect/null-z maps use zero outside transformed valid-center support; p maps '
            'use one. Always use the valid-center mask to distinguish background from '
            'measured zero effects. Binary masks and neighborhood-size counts preserve '
            'native values; T1 voxel counts are not independent tests or native feature counts.\n\n'
            'For display, load `*_delta_z_sig_fwe_p05_T1.nii.gz` on the reference T1. '
            'These effects are masked by the transformed original joint max-stat FWE p < .05 '
            'decisions. Unthresholded effect, null-z, p, support, and binary significance '
            'maps are also included. See projection_manifest.json for hashes and checks.\n')
        publish.rename(destination)
    print(f'Subject {subject}: {len(records)} maps + 6 thresholded overlays verified -> {destination}', flush=True)
    return str(destination)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--subjects', type=int, nargs='+', default=list(range(2, 7)))
    parser.add_argument('--workers', type=int, default=2)
    args = parser.parse_args()
    fsl = Path(os.environ.get('FSLDIR', str(Path.home() / 'fsl')))
    flirt = Path(shutil.which('flirt') or fsl / 'bin' / 'flirt')
    assert flirt.is_file(), f'FSL FLIRT unavailable: {flirt}'
    assert len(set(args.subjects)) == len(args.subjects)
    # Preflight every participant before starting any projections.
    for subject in args.subjects:
        base = args.root / 'MRI' / f'subj_{subject}' / 'nifti'
        assert not (base / 'olf_context_similarity_searchlight/global/native_T1').exists()
        for path in [base / 'anat/betT1brain.nii.gz', base / 'coreg/func2T1mat']:
            assert path.is_file(), path
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        list(executor.map(lambda s: project_subject(s, args.root.resolve(), flirt), args.subjects))


if __name__ == '__main__':
    main()
