# AGENTS.md

This file defines the durable operating rules for agents working in the OX repository. Read it before changing analysis code or interpreting results.

If documents disagree, use this order:

1. `AGENTS.md` for operating rules and analysis guardrails.
2. `project_info.md` for the current study design, canonical inputs/ROIs, completed analyses, and current interpretation.
3. Specialized technical documents for implementation details.
4. `roadmap.md` for mutable planning only; it is not a source of truth.

## Project and repository

OX is a human olfactory fMRI project studying how odor identity and semantic context are represented and combined to generate flexible odor perception.

Repository layout:

- `scripts/`: main runnable MATLAB analysis scripts.
- `utils/`: helper functions and shell utilities.
- `environment/`: environment/setup helpers.
- `ROIs/`: atlas and ROI resources tracked in Git.
- `notes/`: brief reminders and subject-specific notes.
- `archive/`: superseded or discarded code; do not use as the default implementation.
- `behavior/`: behavior result files. including raw and extracted intensity/pleasantness ratings.
- `cuelist/`: raw cue lists ran in actual experiment sessions.
- `data_collection_notes/`: data collection notes collected from actual experiment sessions.
- `labchart/`: raw labchart export files and extracted event files.
- `MRI/`: all nifti files and analysis results directly ran on these files.
- `RDMs/`: derived distance matrices and associated RSA analyses.
- `results/`: some group level results.
- `tmp/`: temporary directory for misc logs and verification records. 

Large raw and derived MRI/data directories are intentionally ignored by Git.  

## Documentation routing

Before specialized work, read the relevant document:

- `project_info.md`: current study design, subject set, canonical inputs, ROI definitions, current analysis status, output locations, and current scientific interpretation.
- `GLMSINGLE_CONFOUNDS.md`: exact construction and provenance rules for run-wise GLMsingle nuisance regressors.
- `NEMO_OX_PIPELINE_COMPARISON.md`: audit of NEMO versus OX and known limitations of legacy OX analyses.
- `data_collection_exclusion_review.md`: acquisition-note review and candidate QC issues. This is not a final exclusion list. Refer to this only when discussing trial exclusion.  
- `roadmap.md`: Deprecated. Do not use. 
- `notes/project_notes.md`: scratch reminders and short subject-specific notes.

## Core study contract

Unless `project_info.md` explicitly records a finalized change:

- Analysis sample is `subj_2` through `subj_6` (`n = 5`).
- Each analyzed participant has 80 acquisition runs, 10 trials per run, and 800 trials total.
- The same 20 odor identities occur in four contexts: `PERSON`, `FOOD`, `LOCATION`, and `CONTROL`.
- There are 200 trials per context.
- `PERSON`, `FOOD`, and `LOCATION` are the three semantic contexts. Do not treat `CONTROL` as a fourth semantic context in analyses explicitly about semantic-context similarity.
- Use the uppercase context labels above in saved metadata/results.
- TR is 0.76 s.

The trial sequence is:

```text
fixation -> context cue -> narrated scenario -> "3, 2, 1, sniff" -> odor/sniff -> ratings
```

Pleasantness and intensity are distinct behavioral variables. 

## Authoritative timing and trial metadata

- `labchart/extracted_events/subjN_events_bm.mat` is authoritative for run ordinal, frame count, MRI onset, sniff TTL provenance, and inhale onset.
- Load exact odor, context, run, and session labels with `OX_load_trial_metadata` rather than rebuilding trial order from filenames or directory order.
- Enumerate functional runs with `OX_discover_functional_runs`; use exact acquisition run IDs as grouping variables.
- Subject 3, runs 1-10 have reversed respiratory polarity. Use `OX_get_respiration_polarity`; do not implement a second ad-hoc correction.
- Never infer run correspondence from filenames when an explicit metadata mapping exists.

For physiology/TTL work, `GLMSINGLE_CONFOUNDS.md` contains the detailed source-of-truth rules. 

## Canonical fMRI inputs

Use the physiology-regressed GLMsingle outputs for new multivariate analyses unless the scientific question explicitly requires another derivative:

- Sniff/odor aligned: `MRI/subj_N/nifti/sniff_single_trial_by_category_physio/TYPED_FITHRF_GLMDENOISE_RR.mat`
- Countdown aligned: `MRI/subj_N/nifti/countdown_single_trial_by_category_physio/TYPED_FITHRF_GLMDENOISE_RR.mat`

Sniff-aligned estimates are primary for odor-perception claims. Countdown-aligned estimates are an earlier comparison/sensitivity analysis. Cue-aligned analyses remain exploratory because the narrated cue is difficult to represent as a single event.

For canonical GLMsingle files:

- `modelmd` rows correspond exactly to `find(gm_mask)`.
- Columns correspond exactly to the trial metadata order.
- Assert both mappings before scoring or reshaping data.

## Confounds and censoring

Before modifying nuisance construction, read `GLMSINGLE_CONFOUNDS.md`.

Non-negotiable conventions:

- Bad fMRI volumes are represented by one-hot spike regressors; volumes are not deleted.
- Run-wise GLMsingle nuisance matrices may contain different numbers of columns because different runs have different numbers of spikes. Do not pad them to equal width.
- Preserve run alignment between motion, physiology, events, and functional images.
- Do not silently alter respiration regressors; respiration is both a scientific variable and an fMRI confound in this project.

## Data safety and provenance

- Treat raw MRI, raw LabChart/respiration/TTL data, manually curated ROI masks, and established preprocessing derivatives as immutable unless explicitly instructed otherwise.
- Do not overwrite an established analysis output when a new output directory can preserve provenance.
- Do not silently drop subjects, runs, trials, frames, odors, or voxels.
- Any exclusion must be explicit and reproducible.
- Preserve subject/session/run/trial identifiers through every processing stage.
- Save enough metadata for a derived result to identify its inputs, subject set, trial set, ROI/restriction criteria, model settings, random seed/permutation count where relevant, and producing script.

`data_collection_exclusion_review.md` is a notes-derived QC review, not an automatic exclusion table. Exclusions would be stated explicitly if needed. 

## ROI and spatial conventions

Use `project_info.md` for the current active ROI sets and analysis-specific restriction masks.

Important guardrails:

- The unqualified phrase **olfactory ROIs** refers to the focused five-ROI set: `PirF`, `PirT`, `AON`, `olfOFC`, `olfAMG`.
- The seven-mask restricted-searchlight support set (`TU`, `AON`, `PirF`, `PirT`, `olfAMG`, `EC`, `HIPP`) is a different object. Spell it out rather than calling it simply “olfactory ROIs.”
- `old_rois/` is an immutable archive for reproducing legacy analyses. Do not use it as the default ROI source for new work.
- Participant-native masks/maps must not be combined voxelwise across subjects before transformation to a common space.
- Report final usable feature counts after all ROI, gray-matter, functional-restriction, and finite-data intersections; anatomical mask size is not the same as final feature count.
- Require at least 10 usable features unless a finalized analysis explicitly defines another threshold. Mark smaller ROIs/neighborhoods as insufficient rather than silently lowering the threshold.

When a new analysis needs functional restriction, follow the current convention recorded in `project_info.md`; at present the preferred future-analysis restriction is the positive Odor > Rest uncorrected `p < .001` mask unless the analysis has a documented reason to use another threshold.

## Multivariate analysis defaults

- Use simple neural distance (Pearson correlation distance, `1 - r`) and linear model fits for RSA and related analyses when applicable.
- Do not rank-transform RDMs unless a specific analysis is explicitly designed to do so.
- Default prediction cross-validation is leave-one-run-out (LORO). Trial-wise leave-one-out is legacy/provisional because same-run information can remain in training.
- Use exact acquisition run IDs for folds; do not infer folds from trial ordinal alone.
- Build templates and fit preprocessing/scaling parameters from training data only when cross-validation requires independence.
- Permutations must rerun the complete scoring/CV pipeline, not only shuffle a final score.
- For decoding, shuffle labels within exact acquisition run so run-level counts and nuisance structure remain fixed.
- Current multivariate analyses center each voxel's GLMsingle betas within run. For semantic-only similarity, compute the centering mean from semantic trials only so `CONTROL` does not enter indirectly. Describe this as transductive within-run normalization.
- Pearson pattern similarity implies spatial demeaning. Current template and split-half analyses do not scale or whiten voxels unless explicitly stated.

Reference chance levels:

- Four-way context classification: 0.25.
- Twenty-way odor classification: 0.05.
- Three-way cross-odor semantic-context classification: 1/3.

Retain participant-level results, raw accuracies/effects, and confusion matrices where relevant even when group summaries use null-centered evidence or accuracy-minus-chance.

## Statistical inference

The dataset is a dense-sampling/high-precision design with five analyzed participants. Population-level random-effects claims require special caution.

- Participant-level effects and plots must remain prominent.
- Do not treat trials, odor pairs, voxels, runs, or repeated split-half realizations as independent substitutes for subjects.
- Avoid conventional one-sample group t-tests as the default inferential strategy for new coefficient-based analyses.
- For new coefficient-based inference, generate within-subject permutation nulls and individual statistics. Construct each group null draw by averaging one null beta from each included subject with equal weights, then compare the observed mean beta with that aggregated null distribution.
- Interpret such group permutation summaries as inference about the measured subjects, not a population random-effects estimate.
- Exact five-subject sign-flip tests contain only 32 sign patterns; do not overstate their resolution.
- Production permutation tests use add-one empirical p-values: `(1 + exceedances) / (N permutations + 1)`.
- State the multiple-comparison family and whether correction uses FDR, a maximum statistic, or another prespecified method.
- Keep smoke tests visibly separate from production results. Small smoke tests are for implementation validation, not inference.

## Legacy analysis warnings

Existing code/results are not automatically approved final analyses.

In particular, earlier whole-brain LOO searchlight context/odor results are provisional because they used trial-wise leave-one-out, same-data R-squared center selection, and fixed-effect Stouffer aggregation. Do not use them for a main scientific claim. See `NEMO_OX_PIPELINE_COMPARISON.md`.

NEMO/Sagar is methodological precedent and a useful comparison, not a pipeline to copy wholesale into OX.

## Coding and change discipline

- Inspect existing runners/utilities before creating a parallel implementation.
- Prefer extending the current canonical implementation over adding near-duplicate code.
- Keep changes scoped to the requested analysis; avoid unrelated refactors during scientific changes.
- Prefer explicit subject/session/run identifiers over filesystem ordering.
- Avoid hard-coded absolute paths when repository setup/configuration already provides the path.
- Do not launch a large recomputation solely to test a small code change when a representative subset or dry run is sufficient.

Before finishing a substantive change, check:

1. Expected subjects, runs, trials, and odors are present.
2. Array dimensions and metadata mappings match.
3. No unintended NaNs, empty ROIs, missing runs, or silent exclusions were introduced.
4. Cross-validation and permutation exchangeability are correct for the scientific question.
5. New outputs are separated from smoke tests and legacy results.
6. Documentation is updated if the change alters a finalized subject set, trial definition, canonical input, ROI set, functional restriction, CV scheme, inference scheme, or output location.
