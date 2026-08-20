#!/usr/bin/env python3
"""Inventory subject-space decoding ROI masks and count their voxels."""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import statistics
import subprocess
from collections import Counter
from pathlib import Path


DEFAULT_SUBJECTS = (2, 3, 4, 5, 6)
EXPECTED_MASKS_PER_SUBJECT = 66


def run_text(*args: str | Path) -> str:
    return subprocess.run(
        [str(arg) for arg in args],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def fsl_voxels(image: Path, mask: Path | None = None) -> tuple[int, float]:
    command: list[str | Path] = ["fslstats", image]
    if mask is not None:
        command.extend(["-k", mask])
    command.append("-V")
    fields = run_text(*command).split()
    if len(fields) != 2:
        raise RuntimeError(f"Unexpected fslstats -V output for {image}: {fields}")
    return int(float(fields[0])), float(fields[1])


def fsl_value(image: Path, key: str) -> float:
    return float(run_text("fslval", image, key))


def check_geometry(image: Path, reference: Path) -> None:
    for key in ("dim1", "dim2", "dim3", "pixdim1", "pixdim2", "pixdim3"):
        image_value = fsl_value(image, key)
        reference_value = fsl_value(reference, key)
        if not math.isclose(image_value, reference_value, rel_tol=0, abs_tol=1e-6):
            raise RuntimeError(
                f"Geometry mismatch for {key}: {image}={image_value}, "
                f"{reference}={reference_value}"
            )


def resolve_manifest_mask(root: Path, roi_dir: Path, manifest_path: str) -> Path:
    saved_path = Path(manifest_path)
    filename = saved_path.name
    if filename.endswith(".nii.gz"):
        alternate_filename = filename.removesuffix(".gz")
    elif filename.endswith(".nii"):
        alternate_filename = f"{filename}.gz"
    else:
        alternate_filename = filename
    candidates = (
        saved_path,
        roi_dir / filename,
        roi_dir / alternate_filename,
    )
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    raise FileNotFoundError(f"ROI mask listed in manifest was not found: {manifest_path}")


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--subjects",
        nargs="+",
        type=int,
        default=list(DEFAULT_SUBJECTS),
        help="Numeric subject IDs (default: 2 3 4 5 6)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=root / "MRI" / "group" / "roi_decoding_inventory",
    )
    args = parser.parse_args()

    for command in ("fslstats", "fslval"):
        if shutil.which(command) is None:
            raise RuntimeError(f"Required FSL command is unavailable: {command}")

    subject_names = [f"subj_{number}" for number in args.subjects]
    records: list[dict[str, object]] = []

    for subject in subject_names:
        coreg_dir = root / "MRI" / subject / "nifti" / "coreg"
        roi_dir = coreg_dir / "roi_decoding"
        manifest = roi_dir / "roi_manifest.tsv"
        gm_mask = coreg_dir / "gm_mask_thr05_func.nii"
        if not manifest.is_file():
            raise FileNotFoundError(f"Missing ROI manifest: {manifest}")
        if not gm_mask.is_file():
            raise FileNotFoundError(f"Missing functional-space GM mask: {gm_mask}")

        subject_records: list[dict[str, object]] = []
        with manifest.open(newline="") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            for manifest_row in reader:
                image = resolve_manifest_mask(root, roi_dir, manifest_row["output_file"])
                check_geometry(image, gm_mask)

                value_range = [float(value) for value in run_text("fslstats", image, "-R").split()]
                if value_range != [0.0, 1.0]:
                    raise RuntimeError(f"ROI mask is not binary: {image} has range {value_range}")

                mask_voxels, mask_volume = fsl_voxels(image)
                gm_voxels, gm_volume = fsl_voxels(image, gm_mask)
                manifest_voxels = int(manifest_row["voxel_count"])
                if mask_voxels != manifest_voxels:
                    raise RuntimeError(
                        f"Manifest mismatch for {image}: saved={manifest_voxels}, "
                        f"measured={mask_voxels}"
                    )

                record: dict[str, object] = {
                    "subject": subject,
                    "roi_id": manifest_row["roi_id"],
                    "hemisphere": manifest_row["hemisphere"],
                    "source": manifest_row["source"],
                    "contributing_labels": manifest_row["contributing_labels"],
                    "mask_file": str(image.relative_to(root)),
                    "mask_voxels": mask_voxels,
                    "gm_overlap_voxels": gm_voxels,
                    "gm_overlap_fraction": round(gm_voxels / mask_voxels, 6),
                    "voxel_volume_mm3": round(mask_volume / mask_voxels, 6),
                    "mask_volume_mm3": round(mask_volume, 3),
                    "gm_overlap_volume_mm3": round(gm_volume, 3),
                }
                subject_records.append(record)

        if len(subject_records) != EXPECTED_MASKS_PER_SUBJECT:
            raise RuntimeError(
                f"Expected {EXPECTED_MASKS_PER_SUBJECT} manifest rows for {subject}; "
                f"found {len(subject_records)}"
            )
        key_counts = Counter(
            (row["roi_id"], row["hemisphere"]) for row in subject_records
        )
        duplicates = [key for key, count in key_counts.items() if count != 1]
        if duplicates:
            raise RuntimeError(f"Duplicate ROI/hemisphere keys for {subject}: {duplicates}")
        records.extend(subject_records)

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    long_fields = [
        "subject",
        "roi_id",
        "hemisphere",
        "source",
        "contributing_labels",
        "mask_file",
        "mask_voxels",
        "gm_overlap_voxels",
        "gm_overlap_fraction",
        "voxel_volume_mm3",
        "mask_volume_mm3",
        "gm_overlap_volume_mm3",
    ]
    write_tsv(output_dir / "decoding_roi_voxel_counts_long.tsv", long_fields, records)

    first_subject = subject_names[0]
    ordered_keys = [
        (str(row["roi_id"]), str(row["hemisphere"]))
        for row in records
        if row["subject"] == first_subject
    ]
    record_by_key = {
        (str(row["roi_id"]), str(row["hemisphere"]), str(row["subject"])): row
        for row in records
    }

    metadata_fields = ["roi_id", "hemisphere", "source", "contributing_labels"]
    for metric, filename in (
        ("mask_voxels", "decoding_roi_voxel_counts_wide.tsv"),
        ("gm_overlap_voxels", "decoding_roi_gm_overlap_counts_wide.tsv"),
    ):
        wide_rows: list[dict[str, object]] = []
        for roi_id, hemisphere in ordered_keys:
            exemplar = record_by_key[(roi_id, hemisphere, first_subject)]
            wide_row: dict[str, object] = {
                field: exemplar[field] for field in metadata_fields
            }
            for subject in subject_names:
                wide_row[subject] = record_by_key[(roi_id, hemisphere, subject)][metric]
            wide_rows.append(wide_row)
        write_tsv(output_dir / filename, metadata_fields + subject_names, wide_rows)

    bilateral_fields = ["roi_id", "source", "contributing_labels"]
    for subject in subject_names:
        bilateral_fields.extend(
            [f"{subject}_mask_voxels", f"{subject}_gm_overlap_voxels"]
        )
    bilateral_rows: list[dict[str, object]] = []
    for roi_id, hemisphere in ordered_keys:
        if hemisphere != "bilateral":
            continue
        exemplar = record_by_key[(roi_id, hemisphere, first_subject)]
        bilateral_row: dict[str, object] = {
            "roi_id": roi_id,
            "source": exemplar["source"],
            "contributing_labels": exemplar["contributing_labels"],
        }
        for subject in subject_names:
            subject_record = record_by_key[(roi_id, hemisphere, subject)]
            bilateral_row[f"{subject}_mask_voxels"] = subject_record["mask_voxels"]
            bilateral_row[f"{subject}_gm_overlap_voxels"] = subject_record[
                "gm_overlap_voxels"
            ]
        bilateral_rows.append(bilateral_row)
    write_tsv(
        output_dir / "decoding_roi_bilateral_voxel_counts.tsv",
        bilateral_fields,
        bilateral_rows,
    )

    summary_rows: list[dict[str, object]] = []
    for roi_id, hemisphere in ordered_keys:
        roi_records = [
            record_by_key[(roi_id, hemisphere, subject)] for subject in subject_names
        ]
        mask_counts = [int(row["mask_voxels"]) for row in roi_records]
        gm_counts = [int(row["gm_overlap_voxels"]) for row in roi_records]
        gm_fractions = [float(row["gm_overlap_fraction"]) for row in roi_records]
        exemplar = roi_records[0]
        summary_rows.append(
            {
                "roi_id": roi_id,
                "hemisphere": hemisphere,
                "source": exemplar["source"],
                "contributing_labels": exemplar["contributing_labels"],
                "n_subjects": len(subject_names),
                "mean_mask_voxels": round(statistics.mean(mask_counts), 1),
                "min_mask_voxels": min(mask_counts),
                "max_mask_voxels": max(mask_counts),
                "mean_gm_overlap_voxels": round(statistics.mean(gm_counts), 1),
                "min_gm_overlap_voxels": min(gm_counts),
                "max_gm_overlap_voxels": max(gm_counts),
                "mean_gm_overlap_fraction": round(statistics.mean(gm_fractions), 3),
            }
        )

    summary_fields = [
        "roi_id",
        "hemisphere",
        "source",
        "contributing_labels",
        "n_subjects",
        "mean_mask_voxels",
        "min_mask_voxels",
        "max_mask_voxels",
        "mean_gm_overlap_voxels",
        "min_gm_overlap_voxels",
        "max_gm_overlap_voxels",
        "mean_gm_overlap_fraction",
    ]
    write_tsv(output_dir / "decoding_roi_voxel_summary.tsv", summary_fields, summary_rows)

    print(f"Validated {len(records)} masks across {len(subject_names)} subjects.")
    print(f"Wrote ROI voxel-count inventory to {output_dir}")


if __name__ == "__main__":
    main()
