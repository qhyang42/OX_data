#!/usr/bin/env python3
"""Run after export_rsa_tsnr_inputs.m; compute descriptive tSNR in saved RSA features."""
import csv
import json
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR', '/tmp/ox_tsnr_matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/ROI_tSNR'
NAMES = ['olf_TU', 'AON', 'PirF', 'PirT', 'olfAMG', 'olfOFC', 'control_GM_putamen', 'control_A1']
LABELS = ['TU', 'AON', 'PirF', 'PirT', 'olfAMG', 'olfOFC', 'Putamen', 'A1']

def write_csv(name, rows):
    with (OUT / name).open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

def main():
    inputs = json.loads((OUT / 'inputs.json').read_text())
    assert [s['subject'] for s in inputs] == list(range(2, 7))
    run_rows, subject_rows, invalid_rows = [], [], []
    values = np.full((5, 8), np.nan)
    for si, sub in enumerate(inputs):
        gm = nib.load(sub['gm_file'])
        rois = sub['rois']
        assert len(sub['runs']) == 80
        indices = np.unique(np.concatenate([r['linear_indices_zero_based'] for r in rois]))
        selections = [np.searchsorted(indices, r['linear_indices_zero_based']) for r in rois]
        voxel_runs = []
        cache = OUT / f"subj_{sub['subject']}_voxel_tsnr.npz"
        if cache.exists():
            saved = np.load(cache)
            np.testing.assert_array_equal(saved['indices'], indices)
            voxel_runs = saved['tsnr']
        else:
            for run in sub['runs']:
                img = nib.load(run['functional_file'], mmap='r')
                assert img.shape[:3] == gm.shape and np.allclose(img.affine, gm.affine, atol=1e-4)
                assert img.shape[3] == run['expected_frames']
                raw = img.dataobj.get_unscaled()
                x = np.asarray(raw.reshape((-1, img.shape[3]), order='F')[indices, :], dtype=np.float64)
                x = x * img.dataobj.slope + img.dataobj.inter
                assert np.isfinite(x).all(), 'Nonfinite time series; investigate source data.'
                sd = x.std(axis=1, ddof=1)
                tsnr = np.divide(x.mean(axis=1), sd, out=np.full(len(indices), np.nan), where=sd > 0)
                voxel_runs.append(tsnr)
                if run['ordinal'] % 20 == 0:
                    print(f"Subject {sub['subject']}: {run['ordinal']}/80 runs", flush=True)
            voxel_runs = np.asarray(voxel_runs)
            np.savez_compressed(cache, indices=indices, tsnr=voxel_runs)
        assert voxel_runs.shape == (80, len(indices))
        for roi, sel in zip(rois, selections):
            data = voxel_runs[:, sel]
            valid = np.isfinite(data).all(axis=0)
            for k, v in zip(*np.where(~np.isfinite(data))):
                invalid_rows.append(dict(subject=sub['subject'], roi=roi['name'],
                    run_id=sub['runs'][k]['id'], linear_index_zero_based=int(indices[sel[v]]),
                    reason='zero temporal SD'))
            n_valid = int(valid.sum())
            j = NAMES.index(roi['name'])
            usable = n_valid >= 10
            for k, run in enumerate(sub['runs']):
                value = float(data[k, valid].mean()) if usable else np.nan
                run_rows.append(dict(subject=sub['subject'], roi=roi['name'], run_id=run['id'],
                    run_ordinal=run['ordinal'], session=run['session'], run_in_session=run['run'],
                    frames=run['expected_frames'], n_rsa_voxels=roi['n_voxels'],
                    n_tsnr_voxels=n_valid, tsnr=value, functional_file=run['functional_file']))
            value = float(data[:, valid].mean(axis=1).mean()) if usable else np.nan
            values[si, j] = value
            subject_rows.append(dict(subject=sub['subject'], roi=roi['name'], n_runs=80,
                n_rsa_voxels=roi['n_voxels'], n_tsnr_voxels=n_valid,
                n_invalid_voxels=int((~valid).sum()), tsnr=value,
                status='usable' if usable else 'insufficient: <10 tSNR features'))
        if sub['subject'] == 5:
            subject_rows.append(dict(subject=5, roi='olf_TU', n_runs=0, n_rsa_voxels=6, n_tsnr_voxels=0, n_invalid_voxels='',
                tsnr='', status='insufficient: 6 < 10 RSA features'))
    assert np.isnan(values[3, 0])
    assert len(run_rows) == 39 * 80
    write_csv('run_tsnr.csv', run_rows)
    write_csv('subject_tsnr.csv', subject_rows)
    if invalid_rows:
        write_csv('invalid_voxel_runs.csv', invalid_rows)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
        'pdf.fonttype': 42, 'svg.fonttype': 'none', 'axes.spines.top': False,
        'axes.spines.right': False})
    fig, ax = plt.subplots(figsize=(11, 5.8))
    bp = ax.boxplot([v[np.isfinite(v)] for v in values.T], positions=np.arange(8),
        widths=.56, patch_artist=True, showfliers=False, whis=1.5,
        medianprops=dict(color='#222222', linewidth=1.6),
        whiskerprops=dict(color='#777777'), capprops=dict(color='#777777'))
    for j, patch in enumerate(bp['boxes']):
        patch.set(facecolor='#d5e4e8' if j < 6 else '#dedede', edgecolor='#727b80')
    colors = ['#0072B2', '#E69F00', '#009E73', '#CC79A7', '#D55E00']
    markers = ['o', 's', '^', 'D', 'v']
    for si in range(5):
        ax.scatter(np.arange(8) + np.linspace(-.17, .17, 5)[si], values[si],
            s=48, color=colors[si], marker=markers[si], edgecolor='white', linewidth=.6,
            zorder=4, label=f'Subject {si+2}')
    ax.axvline(5.5, color='#b9b9b9', linestyle='--', linewidth=1)
    ax.set_xticks(np.arange(8), [f'{label}\n(n={int(np.isfinite(values[:, j]).sum())})' for j,label in enumerate(LABELS)])
    ax.set_ylabel('Temporal signal-to-noise ratio (tSNR)')
    ax.set_ylim(bottom=0)
    ax.set_title('Temporal SNR in RSA regions', loc='left', fontsize=17, pad=48)
    ax.text(0, 1.035, 'Olfactory regions', transform=ax.transAxes, color='#53676d')
    ax.text(.765, 1.035, 'Control regions', transform=ax.transAxes, color='#666666')
    ax.legend(loc='lower left', bbox_to_anchor=(-.005, 1.10), ncol=5, frameon=False, fontsize=9)
    ax.yaxis.grid(True, alpha=.15); ax.set_axisbelow(True)
    fig.text(.09, .045, 'Each point: mean voxelwise tSNR, averaged equally across 80 runs. Boxes summarize subjects.\n'
             'RSA masks + valid tSNR across runs; olfactory: GM + Odor > Rest p < .001; controls: GM only. TU: subject 5 insufficient.',
             fontsize=9, color='#555555')
    fig.subplots_adjust(left=.09, right=.985, bottom=.20, top=.75)
    for ext in ['png', 'pdf', 'svg']:
        fig.savefig(OUT / f'ROI_tSNR.{ext}', dpi=300, facecolor='white')
    plt.close(fig)
    (OUT/'README.md').write_text('''# Temporal SNR in RSA regions

Subjects 2–6; ROIs ordered TU, AON, PirF, PirT, olfAMG, olfOFC, Putamen, A1.
Exact saved `model_feature_indices` are mapped through `find(gm_mask)` into native functional space. These already include finite-beta selection. Olfactory ROIs use GM and positive Odor > Rest uncorrected p < .001; controls use GM only. TU subject 5 is insufficient (6 < 10), matching RSA, and is explicitly unavailable; other cells have five subjects.

For each voxel and acquisition run, tSNR = temporal mean / temporal sample SD (N−1 denominator). Take the arithmetic mean over selected voxels, then average the 80 run means with equal weights per subject. Runs are not concatenated. Inputs are the smoothed, realigned `sr*.nii` time series used as GLMsingle inputs, before GLMsingle physiology/nuisance regression. All frames are retained; no detrending, temporal filtering, nuisance regression, or censoring is added. Thus tSNR includes task-related temporal variance and is preprocessing-dependent; it is not tSNR of GLMsingle beta estimates.

Canonical OX_discover_functional_runs and OX_load_trial_metadata establish run correspondence; frame counts are checked against authoritative extracted events, and spatial shape/affines against native GM. Nonfinite samples cause failure. Zero temporal SD gives undefined tSNR: affected voxels are explicitly listed in invalid_voxel_runs.csv and excluded across all runs of that subject/ROI, producing a fixed valid feature set. At least 10 valid features are required. subject_tsnr.csv reports both original RSA and final tSNR feature counts. No voxelwise combination across subjects.

Boxes: across-subject quartiles, median, whiskers at most 1.5 IQR; all available subject values overlaid with fixed subject offsets/colors/shapes. Descriptive figure, no inferential tests. Functional selection differs between olfactory and control ROIs as in RSA.

Outputs: PNG, vector PDF/SVG; subject_tsnr.csv including final feature counts and insufficient TU entry; run_tsnr.csv with all 3,120 usable subject–ROI–run records; inputs.json with exact voxel indices, source RDM files, masks, and all functional paths. Per-subject voxel tSNR NPZ caches preserve every selected voxel/run, with undefined values represented by NaN. No established analysis outputs overwritten.

Reproduce from repository root: run MATLAB `addpath('scripts'); export_rsa_tsnr_inputs` (requires a new output directory), then Python `scripts/plot_rsa_roi_tsnr.py` with numpy, nibabel, matplotlib. Producing scripts: scripts/export_rsa_tsnr_inputs.m and scripts/plot_rsa_roi_tsnr.py.
''')
    print(values, flush=True)

if __name__ == '__main__':
    main()
