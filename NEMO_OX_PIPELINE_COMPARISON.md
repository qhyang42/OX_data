# NEMO versus OX pipeline audit

**Audit date:** 2026-08-19  
**OX snapshot:** commit `99d26b9914f1b1c1ec7d807542990c905863c332`, plus 28 uncommitted paths  
**NEMO snapshot:** commit `6c4169b7d92b063d9ce011afe8cb4790a4ffe8aa`, clean worktree  
**Scope:** static code review of both repositories, plus read-only checks of OX's existing metadata and result artifacts. No preprocessing or decoding was rerun.

## Bottom line

OX should **not copy NEMO wholesale**. OX's newer ROI decoders are already statistically safer than NEMO's published repository code. The useful lessons from NEMO are upstream: model breathing and artifacts, make the run/session structure explicit, generate QC products, and test whether odor patterns reproduce across genuinely independent sessions.

The current OX whole-brain searchlight maps should be treated as **provisional**. The saved group analysis combines leave-one-trial-out results, so trials from the held-out trial's run remain in training. Searchlight centers are also selected by GLMsingle R² obtained from the same four-category design. In the current group mask, 379/562 context voxels and 557/562 odor-score voxels survive FDR. Such near-ubiquitous odor significance is a red flag for run leakage, circular center selection, residual physiology, or a combination of these—not strong evidence that nearly the entire tested mask contains odor identity.

The four changes most likely to alter the scientific conclusions are:

1. repair and audit the preprocessing execution chain;
2. add respiratory and artifact nuisance regressors before re-estimating trial betas;
3. replace leave-one-trial-out searchlights and same-data R² center selection;
4. perform group inference on subject-level effects, rather than treating Stouffer-combined within-subject p-values as population inference.

## What the two pipelines actually do

| Stage | NEMO | OX | Consequence |
|---|---|---|---|
| Functional alignment | Realigns all runs together to their mean, bridges partial functional coverage to T1 through a whole-brain EPI, then reslices | Uses the same basic bridge, but the runnable script executes batches 1–3 and later 8–9; coregistration/normalization batches 4–6 are not executed top-to-bottom | Existing OX derivatives may be valid if sections were run manually, but the repository cannot reproduce them reliably |
| Smoothing | 2 mm FWHM | 2.5 mm FWHM | Both are light; OX should also retain unsmoothed data for MVPA sensitivity analyses |
| Distortion/slice timing | Neither is implemented | Neither is implemented | Susceptibility distortion is especially important in OFC and piriform cortex |
| Motion nuisance | 24-parameter expansion | 24-parameter expansion | Equivalent baseline motion model |
| Physiology/artifacts | Respiration/spirometry, slice-difference summaries, and one-hot bad-volume regressors | Respiration is loaded for event timing but is not entered into the GLM; no FD/DVARS censoring or tSNR pipeline | OX is more exposed to sniff-locked and motion-linked false patterns |
| Event timing | Behavioral sniff-cue times are adjusted to scanner time per run | LabChart event pulses are aligned to the first MRI pulse, then converted into one global concatenated timeline | OX has the more direct event source, but hard-coded output names and manual exceptions make provenance fragile |
| Response estimation | SPM canonical-HRF or FIR model; outputs voxels × time bins × odors, often averaging repetitions into odor/session patterns | GLMsingle; outputs one beta per trial from a four-category design | OX preserves trial-level information, but trial order and run IDs must be exact |
| Primary pattern analysis | Cross-session odor-pattern similarity and ROI RSA; voxelwise ridge encoding from perceptual descriptors | Template, evidence-score, SVM, nested LDA/SVM, searchlight, and context-transfer decoders | OX asks richer questions and has stronger permutation machinery |
| Cross-validation | Cross-session pattern comparison or leave-one-odor-out encoding; several tuning choices are precomputed or fixed | ROI defaults can be leave-one-run-out; saved group searchlights are leave-one-trial-out | Final OX analyses should be run-grouped everywhere |
| Group inference | Mainly averages three subjects' bootstrap distributions | Stouffer combination of five subjects' voxelwise permutation p-values | Neither implementation is a population random-effects analysis |

## NEMO preprocessing, step by step

### 1. Spatial preprocessing

NEMO builds one SPM batch and runs it in full:

- It realigns whole-brain EPI images and functional runs separately ([make_preprocessing_job.m, lines 37–117](/Users/qhyang/Documents/GitHub/NEMO_scripts/basic_preprocess/make_preprocessing_job.m:37)). Functional runs use register-to-mean (`rtm = 1`).
- It registers the whole-brain EPI to T1, then the partial functional EPI to the whole-brain EPI and applies that transform to the functional images ([lines 120–148](/Users/qhyang/Documents/GitHub/NEMO_scripts/basic_preprocess/make_preprocessing_job.m:120)).
- It estimates T1-to-MNI normalization but leaves functional normalization commented out, preserving native-space analyses ([lines 151–172](/Users/qhyang/Documents/GitHub/NEMO_scripts/basic_preprocess/make_preprocessing_job.m:151)).
- It reslices the functional images and smooths them with a 2 mm FWHM kernel ([lines 174–205](/Users/qhyang/Documents/GitHub/NEMO_scripts/basic_preprocess/make_preprocessing_job.m:174)).

This is the conceptual source of OX's current partial-EPI → whole-brain-EPI → T1 strategy. It is reasonable, but neither repository includes susceptibility-distortion correction, automated registration QC, or a complete machine-readable provenance record.

### 2. Nuisance modeling and timing

NEMO goes well beyond motion-only regression:

- It expands six realignment parameters to 24 columns: motion, derivative, square, and squared derivative ([make_nusiance_regressors.m, lines 35–49](/Users/qhyang/Documents/GitHub/NEMO_scripts/basic_preprocess/make_nusiance_regressors.m:35)).
- It aligns breathing/task recordings to the first scanner trigger and saves a per-run time adjustment ([lines 62–84](/Users/qhyang/Documents/GitHub/NEMO_scripts/basic_preprocess/make_nusiance_regressors.m:62)).
- It smooths, high-pass filters, standardizes, integrates, and downsamples respiratory and spirometry traces to scan resolution ([lines 86–142](/Users/qhyang/Documents/GitHub/NEMO_scripts/basic_preprocess/make_nusiance_regressors.m:86)).
- It computes odd-even slice differences and slice-variance measures, their derivatives/squares, and one-hot regressors for volumes exceeding a within-run z threshold of 5 ([lines 149–229](/Users/qhyang/Documents/GitHub/NEMO_scripts/basic_preprocess/make_nusiance_regressors.m:149)).
- It computes ROI-level temporal SNR as mean signal divided by temporal standard deviation ([t_snr.m, lines 1–88](/Users/qhyang/Documents/GitHub/NEMO_scripts/basic_preprocess/t_snr.m:1)).

The ideas are worth adopting, but the implementation should not be copied literally. It contains manual thresholds and a likely bad-volume indexing error: `bv_counter` accumulates across runs while each run allocates only `max_badvol` columns ([make_nusiance_regressors.m, lines 190–222](/Users/qhyang/Documents/GitHub/NEMO_scripts/basic_preprocess/make_nusiance_regressors.m:190)). OX should implement validated respiratory regressors and FD/DVARS censoring with run-local indexing instead.

### 3. Odor-response construction

NEMO concatenates runs in time, shifts each run's event times by its duration, adds the concatenated nuisance matrix, and then calls `spm_fmri_concatenate` to restore run boundaries ([FIR_model.m, lines 117–184](/Users/qhyang/Documents/GitHub/NEMO_scripts/FIR_model.m:117)). Each odor/session becomes a condition. The basis is either a canonical HRF or an FIR basis ([lines 143–170](/Users/qhyang/Documents/GitHub/NEMO_scripts/FIR_model.m:143)). The beta maps are finally packed into a `voxel × basis/time × odor` tensor ([lines 203–263](/Users/qhyang/Documents/GitHub/NEMO_scripts/FIR_model.m:203)).

The important design principle is separation of:

1. raw time-series preprocessing;
2. nuisance-aware response estimation;
3. a compact, explicitly indexed pattern tensor.

OX currently reaches a similar endpoint through GLMsingle, but the ordering contract is implicit: the decoder assumes GLMsingle column *i* matches cue-list label *i* and often infers 80 consecutive runs of 10 trials.

## How NEMO constructs its “decoders”

NEMO does not contain one unified classifier. It contains three related pattern-analysis constructions.

### A. Cross-session odor identity (“Mnemonic”)

For each anatomical ROI, NEMO:

1. intersects the ROI with gray matter and an odor-responsive functional mask;
2. loads independent session-wise odor patterns;
3. correlates every session-1 odor pattern with every session-2 odor pattern across voxels;
4. treats the matrix diagonal as same-odor similarity and off-diagonal cells as different-odor similarity;
5. bootstraps odors and averages across three session pairs and three subjects.

The construction is visible at [mnemonic.m, lines 38–117](/Users/qhyang/Documents/GitHub/NEMO_scripts/mnemonic.m:38). A small helper can turn a square similarity matrix into top-1 identity accuracy by asking whether each diagonal element is the column maximum ([decoder_mat.m, lines 1–17](/Users/qhyang/Documents/GitHub/NEMO_scripts/common_functions/decoder_mat.m:1)), although the main mnemonic script reports similarity contrast rather than this classification accuracy.

This is a useful validation idea for OX because the train and test patterns come from independent sessions. However, NEMO's reported formula appears inconsistent: it Fisher-transforms the diagonal with `atanh` but transforms the off-diagonal with `tanh` ([mnemonic.m, line 113](/Users/qhyang/Documents/GitHub/NEMO_scripts/mnemonic.m:113)). Its “permutations” in that script are odor bootstraps, not a null label permutation. These details should not be copied.

### B. Perceptual encoding model

NEMO's EEM is an **encoding model**, not a class decoder. For each voxel, it predicts odor response from behavioral perceptual descriptors:

1. choose one FIR time bin per voxel from a precomputed `argmaxes2.mat`;
2. leave one odor out;
3. standardize neural responses using the training odors only;
4. fit voxelwise ridge regression from descriptors to response;
5. predict the held-out odor response;
6. correlate predicted and observed response across held-out odors to obtain a voxel map.

See [EEM_LOOCV.m, lines 96–157](/Users/qhyang/Documents/GitHub/NEMO_scripts/EEM_LOOCV.m:96). Training-only standardization is good. The default fixed lambda and precomputed time-bin selection are less safe; if either was chosen using all odors, performance is optimistic. OX's nested decoder already handles hyperparameter selection more rigorously and should keep that approach.

### C. ROI RSA

NEMO also correlates neural odor-similarity matrices with perceptual and chemical similarity matrices, optionally regressing task/run structure and bootstrapping odor pairs ([RSA_ROI.m, lines 116–238](/Users/qhyang/Documents/GitHub/NEMO_scripts/RSA_ROI.m:116)). This is relevant to OX's planned behavior link, but it is complementary to decoding, not evidence of classification.

## What OX already does better

The newer OX decoder cores should be preserved:

- ROI masks are checked for matching dimensions and affine geometry, and GLMsingle rows are asserted to match gray-matter-mask indices ([OX_roi_decode_core.m, lines 175–229](/Users/qhyang/Desktop/OX_DATA/utils/OX_utilities/OX_roi_decode_core.m:175)).
- Template, evidence, and SVM models are explicit; training-fold scaling is applied to the held-out fold without using its statistics ([OX_roi_decode_core.m, lines 645–699](/Users/qhyang/Desktop/OX_DATA/utils/OX_utilities/OX_roi_decode_core.m:645)).
- Labels are permuted within runs and the entire preprocessing/CV pipeline is rerun; ROI-wise FDR and max-statistic FWE are reported ([OX_roi_decode_core.m, lines 296–360](/Users/qhyang/Desktop/OX_DATA/utils/OX_utilities/OX_roi_decode_core.m:296)).
- The focused olfactory decoder nests model/hyperparameter choice inside grouped outer folds and repeats nested selection under permutation ([OX_olfroi_decode_core.m, lines 220–306](/Users/qhyang/Desktop/OX_DATA/utils/OX_utilities/OX_olfroi_decode_core.m:220)).
- Cross-context and cross-odor transfer analyses address invariance more directly than NEMO's general identity comparison.

These are strong foundations. The main risks enter before the decoder or through the current runner settings.

## Prioritized improvements for OX

### P0 — required before interpreting the main maps

#### 1. Make preprocessing one auditable, top-to-bottom function

**Problem.** `fmri_preproc.m` runs batches 1–3 ([lines 28–94](/Users/qhyang/Desktop/OX_DATA/scripts/fmri_preproc.m:28)), defines coregistration in slots 4–5, overwrites those same slots with reverse-direction registrations ([lines 96–167](/Users/qhyang/Desktop/OX_DATA/scripts/fmri_preproc.m:96)), defines normalization in slot 6, and finally runs only slots 8–9 ([lines 170–250](/Users/qhyang/Desktop/OX_DATA/scripts/fmri_preproc.m:170)). Therefore a full script run does not execute coregistration or normalization. Manual section execution may have produced the current files, but that cannot be inferred from the code.

**Action.** Replace the script with `OX_preprocess_subject(config, subjidx)` that:

- resolves and sorts every run explicitly;
- runs one unambiguous transform direction for whole-brain EPI → T1 and partial EPI → whole-brain EPI;
- saves the exact input list, SPM batch, transform files, software versions, and completion status;
- stops if any expected run, mean image, affine, or output is missing;
- writes registration overlay PNGs for partial EPI/whole-brain EPI/T1 and ROI/functional intersections.

**Acceptance test.** A clean derivative directory can be generated with one command; rerunning it produces the same file manifest; every subject passes affine/dimension checks and human review of the overlays.

#### 2. Rebuild nuisance regression around olfactory physiology

**Problem.** OX's nuisance file contains exactly 24 motion columns for subjects 2–6. This comes from motion, derivative, square, and squared derivative only ([make_motion_regressors.m, lines 19–35](/Users/qhyang/Desktop/OX_DATA/utils/OX_utilities/make_motion_regressors.m:19)). Respiration is extracted in LabChart preprocessing but never reaches GLMsingle. Sniff timing and respiratory depth can covary with odor identity, context, attention, and head motion.

**Action.** Generate a run-wise confound table containing:

- the existing 24 motion terms;
- framewise displacement and DVARS;
- one-hot spike regressors for volumes crossing preregistered thresholds;
- respiratory phase terms and a slow respiratory-volume/amplitude regressor, plus derivatives where justified;
- non-steady-state indicators and an explicit list of excluded runs/trials.

Pass those run-wise matrices through `opt.extraregressors` when estimating GLMsingle betas. Do not copy NEMO's fixed `z > 5` slice rule without validation. Compare the current and expanded models using tSNR, residual DVARS, beta reliability, and the decoder's null distribution.

**Acceptance test.** Every confound matrix has exactly one row per acquired volume, no NaN/Inf columns, no constant columns after run-wise construction, and a QC report showing censored volumes and respiration/condition correlations.

#### 3. Rerun searchlights with grouped CV and non-circular center masks

**Problem A: folds.** `searchlight_simple.m` requests leave-one-out ([lines 7–17 and 32–43](/Users/qhyang/Desktop/OX_DATA/scripts/searchlight_simple.m:7)). The saved group manifest confirms that all five subjects' current context and odor maps come from `_loo` directories. Trial-wise LOO leaves nine same-run trials in the training set, so run/session fingerprints and temporally correlated estimation noise can support prediction.

**Problem B: center selection.** Searchlight centers are always `R2 > R2Threshold`, even when `RestrictFeaturesToR2=false` ([OX_searchlight_decode_core.m, lines 184–213](/Users/qhyang/Desktop/OX_DATA/utils/OX_utilities/OX_searchlight_decode_core.m:184)). Here R² comes from GLMsingle's four-category design, so the same data and category structure influence which locations are tested. Label permutations do not recompute GLMsingle or R², so they do not reproduce this selection step.

**Action.** Split `CenterMask` from `FeatureMask` and make the final defaults:

- `CrossValidation = 'leave-one-run-out'` or grouped 10-fold with whole runs assigned to folds;
- centers = all valid gray-matter voxels within a prespecified coverage mask, not same-data R²;
- features = all finite gray-matter neighbors;
- optional independent localizer only if its contrast is orthogonal to context/odor identity;
- exact run IDs from `OX_load_trial_metadata`, not inferred equal blocks.

Run the current LOO/R² analysis only as a labeled sensitivity comparison.

**Acceptance test.** A null simulation and shuffled-label run center near chance; LORO effects are reported beside LOO effects; the searchlight coverage no longer depends on category-model R².

#### 4. Change the group claim and group test

**Problem.** `group_searchlight_stouffer.m` combines within-subject p-values with unweighted Stouffer Z ([lines 1–58](/Users/qhyang/Desktop/OX_DATA/scripts/group_searchlight_stouffer.m:1)). This tests combined evidence in these five datasets; it does not estimate between-subject variation and is not a population random-effects test.

**Action.** Use each subject's accuracy-minus-chance or evidence-minus-null map as the group input. Run a one-sample sign-flip/permutation test across subjects, with max-statistic or cluster correction defined at the group level. With only five usable subjects, report individual effects and confidence intervals prominently and label inference exploratory; only 32 sign patterns exist, so exact group resolution is intrinsically limited.

**Acceptance test.** The group result is reproducible from five subject effect maps, reports effect size and subject dispersion, and does not use warped subject p-values as the population-level outcome.

### P1 — high-value robustness work

#### 5. Make one canonical trial manifest

`OX_load_trial_metadata.m` already creates exact session, run, trial, odor, and context fields and validates 800 trials/80 runs ([lines 43–114](/Users/qhyang/Desktop/OX_DATA/utils/OX_utilities/OX_load_trial_metadata.m:43)). Use it everywhere. Retire direct label loading through the absolute path in `OX_get_odor.m` ([lines 1–17](/Users/qhyang/Desktop/OX_DATA/utils/OX_utilities/OX_get_odor.m:1)). Join LabChart events, cue lists, NIfTI runs, motion rows, confounds, and GLMsingle beta columns by explicit IDs.

The present LabChart script also saves hard-coded `subj4` filenames regardless of `subjidx` ([labchart_preproc.m, lines 126–127 and 195–198](/Users/qhyang/Desktop/OX_DATA/scripts/labchart_preproc.m:126)). Replace these with `subjname`, remove silent `try/catch` skipping, and write an exception table rather than inline subject-specific edits.

**Acceptance test.** For every beta column, one manifest row identifies subject/session/run/trial/event time/odor/context; joins are one-to-one; all volumes, trials, and labels reconcile before GLMsingle starts.

#### 6. Address OFC/piriform distortion and registration quality

Neither repository performs susceptibility-distortion correction. For new acquisitions, collect reverse-phase-encoding or field-map data and use them. For existing data, test an appropriate fieldmap-less method, but do not assume it can recover dropout. In all cases, report:

- EPI–T1 boundary overlays around OFC, piriform cortex, amygdala, and temporal pole;
- ROI coverage, gray-matter overlap, mean EPI signal, tSNR, and dropout fraction per subject;
- whether each ROI passes a prespecified coverage threshold.

The current inventory shows some small/variable masks—for example bilateral olfactory tubercle has as few as 15 gray-matter-overlap voxels and bilateral PAC as few as 27. These require stability plots and should not be interpreted from a single threshold.

#### 7. Add independent-session reliability as a companion analysis

Borrow NEMO's best decoder idea: aggregate independent run/session halves, build odor or context templates in one half, and test them in the other. Repeat across balanced split-halves and report reliability alongside LORO accuracy. This asks whether the spatial code reproduces across acquisition blocks, rather than merely whether a trial resembles a template containing trials acquired minutes away.

Use proper Fisher transformation for all correlations and obtain the null by shuffling labels within the appropriate run/session blocks. Do not use NEMO's `atanh(diagonal) - tanh(off-diagonal)` expression.

#### 8. Predeclare decoder and preprocessing sensitivity analyses

Use one primary model and keep alternatives diagnostic:

- Primary: LORO correlation-template or regularized linear model, balanced accuracy, no same-data R² feature selection.
- Sensitivity: unsmoothed versus 2.5 mm smoothed betas; run-centering on/off; training-fold voxel scaling on/off; equalized ROI voxel counts; template versus nested LDA/SVM.
- Generalization: context decoder across odors and odor decoder across contexts, which are closer to the scientific claims than ordinary within-distribution classification.

OX currently subtracts each run's voxel mean using the entire run before CV ([OX_roi_decode_core.m, lines 140–156 and 550–555](/Users/qhyang/Desktop/OX_DATA/utils/OX_utilities/OX_roi_decode_core.m:140)). Because this also uses the held-out run's unlabeled trials, describe it as transductive run normalization and include an analysis without it. The permutation remains internally matched, but the reported accuracy answers a different question from prediction of a single future trial.

### P2 — reproducibility and maintainability

- Move subject IDs, run counts, paths, TR, smoothing kernel, event choice, and thresholds into one versioned config.
- Convert section-driven scripts into functions with inputs/outputs and unique derivative directories.
- Save git commit, dirty-state diff hash, MATLAB/SPM/GLMsingle versions, RNG seed, exact options, input checksums, and trial manifest with every result.
- Add small synthetic tests for event/run alignment, fold isolation, template scoring, permutation exchangeability, and mask-to-beta indexing.
- Expand the README from its current directory list into a runnable dependency graph with one command per stage.

## Recommended execution order

Do not begin by tuning the classifier. Execute the work in this order:

1. **Freeze current results.** Mark existing `_loo` searchlights and their Stouffer group maps as legacy/provisional; do not overwrite them.
2. **Create the canonical manifest.** Reconcile 800 trials, 80 runs, NIfTI volumes, event times, cue-list labels, and motion rows for subjects 2–6. Subject 1 currently has no GLMsingle result and should be explicitly excluded with a reason.
3. **Audit spatial transforms.** Reconstruct which transforms produced every functional mask and inspect EPI/T1/ROI overlays. Re-preprocess only after the transform chain is unambiguous.
4. **Create expanded confounds and QC.** Add respiration and censoring, then regenerate GLMsingle betas in a new derivative directory.
5. **Run decoder smoke tests.** Use two ROIs, zero or a few permutations, serial execution, and exact run labels. Confirm chance behavior with shuffled labels.
6. **Run final subject analyses.** LORO ROI and searchlight models; within-run full-pipeline permutations; independent-session reliability; context/odor transfer.
7. **Run group analysis.** Use subject effect maps/tables and a subject-level random-effects or sign-flip framework. Keep individual-subject plots visible.
8. **Compare old versus new.** Attribute changes separately to grouped CV, center-mask choice, physiological denoising, and spatial preprocessing.

## Concrete run template after the P0 edits

The following is the intended interface for the refactored pipeline; these entry points do not exist yet and should be implemented as part of P0:

```matlab
addpath('/Users/qhyang/Desktop/OX_DATA/environment');
setup_ox;

cfg = OX_load_config('config/analysis.yml');

for subjidx = 2:6
    manifest = OX_build_subject_manifest(cfg, subjidx);
    OX_validate_subject_manifest(manifest);
    OX_preprocess_subject(cfg, subjidx);
    OX_build_confounds(cfg, subjidx);
    OX_estimate_single_trials(cfg, subjidx, 'Alignment', 'odor');

    run_ids = manifest.run_id;
    OX_roi_decode_core(subjidx, 'context', ...
        'CrossValidation', 'leave-one-run-out', ...
        'RunLabels', run_ids, ...
        'RestrictFeaturesToR2', false, ...
        'NumPermutations', 1000);

    OX_searchlight_context(subjidx, ...
        'CrossValidation', 'leave-one-run-out', ...
        'RunLabels', run_ids, ...
        'CenterMask', 'prespecified_gray_matter_coverage', ...
        'RestrictFeaturesToR2', false, ...
        'NumPermutations', 1000);
end

OX_group_effect_permutation(cfg, 2:6, 'context');
```

Before the long run, use `NumPermutations = 2`, two ROIs or a small center subset, and `UseParallel = false`. The final job should refuse smoke-test options and should write a manifest proving the final settings.

## Decision guide

| Question | Recommendation now |
|---|---|
| Can the existing ROI LORO results be explored? | Yes, cautiously. Their decoder implementation is strong, but the betas still lack physiological/artifact regressors and depend on an unreproducible preprocessing chain. |
| Can the existing group searchlight maps support a main claim? | No. They are LOO, use same-data R²-selected centers, and use fixed-effect Stouffer aggregation. |
| Should OX adopt NEMO's FIR estimator? | Not by default. GLMsingle is appropriate for trial-level decoding. Add FIR/canonical-HRF estimates as a timing/reliability sensitivity analysis if the scientific question needs them. |
| Should OX adopt NEMO's decoder code? | No. Keep OX's newer grouped/nested/permutation decoders. Borrow the independent-session comparison concept only. |
| What should be done first? | Fix provenance/manifest/preprocessing execution, then re-estimate betas with respiratory and artifact confounds, then rerun grouped decoders. |

## Audit evidence from current OX artifacts

Read-only checks found:

- subjects 2–6 each have 800-trial odor-aligned and countdown-aligned GLMsingle outputs;
- each has 80 cue-list runs of 10 trials;
- motion files match the total number of functional volumes and contain 24 columns;
- subject 1 lacks a GLMsingle result;
- the current group searchlight manifest explicitly points to `searchlight_context_loo` and `searchlight_score_odor_loo` inputs;
- the group-valid mask contains only 562 voxels, of which 379 context and 557 odor-score voxels are FDR-significant under the current Stouffer procedure.

These checks establish internal dimensions, not correctness of event-to-beta alignment or preprocessing provenance. Those require the manifest and transform audit above.
