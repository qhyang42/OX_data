# OX_DATA

Working repository for the OX human olfactory fMRI project.

## Repository layout

- `scripts/`: main runnable MATLAB analysis scripts
- `utils/`: helper functions and shell utilities
- `environment/`: environment and setup helpers
- `ROIs/`: atlas and ROI resources tracked in Git
- `notes/`: brief project notes and reminders
- `archive/`: discarded or superseded scripts; use only for deliberate legacy reproduction

Large raw and derived data directories are intentionally ignored by Git.

## Documentation map

The Markdown files have distinct roles. Keep those roles separate so current methods, historical audits, and planning do not become mixed together.

### `AGENTS.md` — operating rules

Read this first when using Codex or another coding agent. It contains the durable data contracts, analysis guardrails, statistical defaults, provenance rules, and pointers to specialized documentation.

Update `AGENTS.md` only when a rule should apply broadly to future work in the repository.

### `project_info.md` — current project state

This is the detailed source of truth for the current OX study and analysis state. It records the design, canonical subject-level inputs, ROI definitions, current analysis-specific conventions, completed analyses, output locations, and current interpretation.

Update it when a finalized analysis changes the subject set, trial definition, canonical inputs, ROI set, restriction mask, cross-validation/inference scheme, or output location. Current results and analysis-specific details belong here rather than in `AGENTS.md`.

### `GLMSINGLE_CONFOUNDS.md` — nuisance-regressor implementation

Technical specification for the run-wise GLMsingle confound builder: authoritative inputs, respiratory processing, TTL/window rules, bad-volume handling, final matrix construction, and compatibility details.

Keep implementation detail here; only durable non-negotiable rules should be duplicated in `AGENTS.md`.

### `NEMO_OX_PIPELINE_COMPARISON.md` — historical/methodological audit

Snapshot audit comparing the NEMO and OX pipelines and documenting limitations of older OX analyses. Treat it as methodological history and rationale, not as the current pipeline specification.

### `data_collection_exclusion_review.md` — acquisition-note QC review

Review of handwritten acquisition notes and candidate data-quality issues. It is not a final exclusion list. Final handling must be confirmed against authoritative behavioral, imaging, LabChart/respiration, motion, or acquisition records.

### `roadmap.md` — mutable planning

Figure-oriented analysis planning and next-step notes. This file may contain abandoned, provisional, or superseded ideas. It is not a source of truth for approved methods or current results.

### `notes/project_notes.md` — scratch notes

Short reminders and subject-specific notes that do not belong in the formal project documentation.

## Documentation rule of thumb

Use `AGENTS.md` for **how work should be done**, `project_info.md` for **what the current project/analysis actually is**, specialized Markdown files for **technical or historical detail**, and `roadmap.md` for **what may be done next**.
