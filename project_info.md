# OX experiment and current analysis state

This document records the current OX fMRI study design, canonical analysis inputs and ROIs, completed analyses, current interpretation, and output locations. Durable operating rules and analysis guardrails live in `AGENTS.md`. Specialized implementation details remain in their dedicated technical documents. Paths are relative to the project root unless an absolute path is shown.

## Scientific question

How are odor identity and semantic context represented and combined to generate flexible odor perception in the human brain?

The working figure-level questions are:

1. Where are odor- and context-evoked responses?
2. Where can odor identity and context be decoded during odor perception?
3. Does semantic context change the representation of odor identity?
4. Do neural changes track context-related changes in pleasantness?

## Design, participants, and labels

- Completed analysis sample: `subj_2` through `subj_6` (`n = 5`). `subj_1` is not part of the completed fMRI analysis sample and does not have the canonical GLMsingle derivatives.
- Each completed participant has 80 acquisition runs, 10 trials per run, and 800 trials total.
- The same 20 odor identities occur in four contexts: `PERSON`, `FOOD`, `LOCATION`, and `CONTROL` (no semantic scenario).
- There are 200 trials per context. Analyses restricted to the three semantic contexts therefore use 600 trials.

A trial consists of a context cue and narrated scenario, followed by the auditory countdown/sniff prompt, odor sampling, and pleasantness and intensity ratings:

```text
fixation -> context cue -> narrated scenario -> "3, 2, 1, sniff" -> odor/sniff -> ratings
```

The primary behavioral result to date is that semantic context modulates odor pleasantness. Intensity is retained as a separate measure and should not be substituted for pleasantness.

## Canonical subject-level inputs

### Timing and trial metadata

- BreathMetrics-derived event files are in `labchart/extracted_events/subjN_events_bm.mat`; they contain the canonical run ordinal, frame count, MRI onset, sniff TTL provenance, and inhale onset.
- Exact odor, context, run, and session labels are loaded with `OX_load_trial_metadata`.
- Functional runs are enumerated with `OX_discover_functional_runs`, which sorts numeric session/run IDs and checks the expected run count.
- TR is 0.76 s.
- Subject 3, runs 1-10 have reversed respiratory polarity; the shared correction is in `OX_get_respiration_polarity`.

### Preprocessing and nuisance regression

- Functional preprocessing and native functional-space registration are complete for subjects 2-6.
- Canonical run-wise confounds are in `MRI/subj_N/nifti/glmsingle_confounds/`.
- Each confound matrix contains 24 motion terms, three within-run standardized respiratory terms (`airflow`, squared airflow, and integrated respiratory volume), and one-hot spike regressors for volumes flagged by FD, robust DVARS, or incomplete final-volume respiration.
- `GLMSINGLE_CONFOUNDS.md` records the full construction and provenance rules.

### Canonical single-trial estimates

The canonical physiology-regressed GLMsingle derivatives are:

- Odor/sniff aligned: `MRI/subj_N/nifti/sniff_single_trial_by_category_physio/TYPED_FITHRF_GLMDENOISE_RR.mat`
- Countdown aligned: `MRI/subj_N/nifti/countdown_single_trial_by_category_physio/TYPED_FITHRF_GLMDENOISE_RR.mat`

Both exist for subjects 2-6 and contain 800 trial estimates. The sniff-aligned estimates are primary for odor-perception claims. Countdown-aligned estimates are a temporally earlier comparison/sensitivity analysis. Cue-aligned analyses are exploratory because the extended narrated cue is difficult to model as a single event.

In these files, `modelmd` rows correspond exactly to `find(gm_mask)` and columns correspond to trials in the exact metadata order.

## Mask and ROI conventions

### Gray matter and functional response masks

The subject-native functional-space gray-matter mask is:

```text
MRI/subj_N/nifti/coreg/gm_mask_thr05_func.nii
```

The two saved Odor > Rest restriction masks are:

```text
MRI/subj_N/nifti/first_level_model_sniff_physio/spmT_0001_uncorrected_p001.nii
MRI/subj_N/nifti/first_level_model_sniff_physio/spmT_0001_FWE_p001.nii
```

These come from the physio-regressed sniff-only SPM model. They are thresholded versions of the positive `Odor > implicit rest` contrast, not the generic SPM `mask.nii` estimation mask.

The current analysis-specific convention is:

| Analysis family | Functional restriction |
| --- | --- |
| Standard ROI context/odor decoding (`OX_roi_decode_core`) | `spmT_0001_FWE_p001.nii` |
| Odor-context and cross-odor context-transfer ROI analyses | `spmT_0001_FWE_p001.nii` |
| ROI semantic-context split-half similarity | `spmT_0001_uncorrected_p001.nii` |
| Restricted olfactory-network searchlight | `spmT_0001_uncorrected_p001.nii` |

Final ROI-analysis features are normally:

```text
anatomical ROI AND gm_mask_thr05_func AND selected Odor>Rest mask AND finite GLMsingle rows
```

The audit table for anatomical, GM-overlap, and final feature counts under both functional thresholds is `MRI/group/roi_decoding_inventory/roi_functional_restriction_voxel_counts.csv`.

### Focused olfactory ROIs

The focused five-ROI set used by the ROI semantic-context similarity analysis is:

```text
PirF, PirT, AON, olfOFC, olfAMG
```

- `PirF`: frontal piriform cortex.
- `PirT`: temporal piriform cortex.
- `AON`: anterior olfactory nucleus.
- `olfOFC`: the small olfactory OFC mask, not `fullOFC` or Brainnetome `bn_ofc`.
- `olfAMG`: the union of `MeA + ACo + CeA + PAC`, not whole amygdala and not `nonolfAMG`.

The restricted context-similarity searchlight uses a different, explicitly named extended olfactory-network support set:

```text
TU, AON, PirF, PirT, olfAMG, EC, HIPP
```

This seven-mask union defines both eligible features and searchlight centers after intersection with GM, the uncorrected `p < .001` Odor > Rest mask, and finite GLMsingle rows. It includes olfactory tubercle (`TU`), entorhinal cortex (`EC`), and hippocampus (`HIPP`) and does not include `olfOFC`.

### Active ROI directories

All active decoding masks are bilateral, binary, in each participant's native functional space, and use the suffix `_bilateral_func_thr02.nii`:

```text
MRI/subj_N/nifti/coreg/roi_decoding/primary/
MRI/subj_N/nifti/coreg/roi_decoding/secondary/
```

The active primary set has 13 masks:

```text
AON, PirF, PirT, olfAMG, olfOFC, fullOFC, EC,
maPFC, HIPP, paraHIPP, TP, insula, nACC
```

The active secondary set has 6 masks:

```text
bn_dlpfc, bn_lateral_frontopolar_pfc, bn_ofc,
bn_amygdala, olf_TU, nonolfAMG
```

`nonolfAMG` is `BMA + BLA + LA + PCo`. `old_rois/` is the archive used to reproduce legacy analyses. Each active directory contains `roi_manifest.tsv` with source atlas, contributing labels, hemisphere, output file, and voxel count. Current selection options are `ROISelection='primary'`, `'secondary'`, or `'all'`; `'old'` reproduces archived results.

ROI registration uses the subject transform chain `standard -> T1 -> whole brain -> mean functional`, trilinear interpolation, threshold 0.2, and binarization. Geometry must match the subject's functional reference exactly.

Durable multivariate, cross-validation, permutation, feature-count, and group-inference rules are defined in `AGENTS.md` rather than duplicated here.

## Completed analyses and current interpretation

### Behavioral and univariate work

- Behavioral preprocessing and analysis are complete for the five analyzed participants; context modulation of pleasantness is established at the descriptive/project level.
- `scripts/behavior_analysis.m` regenerates both legacy correlation boxplot variants for pleasantness and intensity for subjects 2–6 in `results/behavior/`, with trial metadata, missing-rating counts, condition-level descriptive statistics, and reproducible sampled correlations. These shared-trial resampling distributions are descriptive, not independent observations for inference.
- Sniff/odor-aligned simple response maps and context contrasts have been generated.
- Cue-aligned context maps were explored and were difficult to interpret. Countdown-aligned maps were more plausible and qualitatively resembled sniff-aligned maps.
- Physio-regressed sniff-only and sniff-by-category first-level SPM models are available for subjects 2-6. The sniff-only model supplies the current functional restriction masks.
- GLMsingle context contrast maps for sniff and countdown estimates are implemented by `scripts/run_glmsingle_context_contrasts.m`.

### Group univariate contrast displays in MNI space

- `scripts/group_fwe_contrasts_mni.py` generates descriptive group displays in `results/contrast_map/` from the existing sniff-aligned, physiology-regressed positive SPM FWE `p < .001` maps for subjects 2–6: Odor > Rest and each of PERSON, FOOD, LOCATION > mean of the other three contexts (including CONTROL).
- Existing `inverse(std2T1mat) × func2T1mat` affine transforms project maps to the repository MNI152 1-mm reference with nearest-neighbour interpolation. Mean thresholded t maps use a fixed denominator of five and zero background; integer significance-count maps and fractions retain individual threshold decisions. SPM estimation-mask coverage counts and individual MNI outputs are also saved.
- These are descriptive averages and overlap maps, not group-level FWE inference. The output README and manifest record source hashes, transforms, validation, coverage, and color conventions.

### ROI organization and QC

- Atlas and in-house masks have been registered to native functional space and reorganized into the active 13-mask primary and 6-mask secondary bilateral sets.
- ROI manifests and group voxel-count inventories have been generated. Geometry checks passed for the focused ROI restriction inventory.

### ROI context and odor decoding

- LORO nearest-template context decoding is complete for the primary ROI set using both sniff- and countdown-aligned physiology-regressed estimates, 1,000 within-run permutations, and participant-level inference.
- LORO 20-way odor decoding is complete for the primary ROI set using countdown-aligned physiology-regressed estimates.
- Plotting code for participant accuracies and context recall is in `scripts/plot_ROI_decoding_results.m`.
- Earlier whole-brain LOO searchlight context/odor maps are provisional; they used trial-wise LOO, same-data R-squared center selection, and fixed-effect Stouffer aggregation. `NEMO_OX_PIPELINE_COMPARISON.md` records the audit.

### Context modulation and generalization analyses

- Odor-context template transfer and cross-odor semantic-context transfer are implemented with LORO training, within-run permutations, synthetic tests, and exact group sign-flip summaries. The currently saved group outputs used `ROISelection='old'` and should be treated as exploratory/legacy until rerun with an explicitly selected active ROI set.
- The primary odor-context transfer result asks whether odor identity is represented more similarly within the same semantic context than across semantic contexts, and also tests transfer from `CONTROL` templates.
- The cross-odor analysis asks whether context patterns learned from 19 odors classify semantic context for the held-out odor.

### Semantic-context split-half similarity

- The focused ROI analysis is complete for subjects 2-6 using the five olfactory ROIs, sniff-aligned physiology-regressed betas, 200 session-balanced run splits, and 5,000 within-run permutations.
- Each split requires all 20 odors in each of `PERSON`, `FOOD`, and `LOCATION` in both halves. Templates average repetitions within odor first and then weight the 20 odors equally.
- The statistic is Fisher-z same-context cross-half pattern similarity minus the mean of the four bidirectional different-context cells.
- The clearest current group-level pattern is stable context-specific representation in `olfOFC`; the group mean effect is also positive in the other four focused ROIs, with more participant and context variability.
- Results are in `MRI/group/context_split_half_similarity_physio/`.

### Formal omnibus and rating-displacement RSA

- TU-only omnibus RSA is saved in `RDMs/omnibus_RSA_TU/`, with new TU RDMs in `RDMs/neural_TU/` and entry point `scripts/run_omnibus_rsa_TU.m`. It uses subjects 2, 3, 4, and 6; subject 5 has 6 usable voxels, below the unchanged 10-voxel minimum. The 5,000 permutations match the original mappings for included subjects. BH-FDR is across the four TU predictors; no predictor passes q < .05 (context beta = -0.0211, p = .0252, q = .1008). Existing ROI inputs/results were preserved.
- Entry points are `scripts/run_omnibus_rsa.m` and `scripts/run_rating_displacement_rsa.m`; shared implementation is in `utils/OX_utilities/`. Both reuse saved neural and behavioral RDMs for subjects 2–6 and the five olfactory ROIs, including all four contexts, with simple neural distance and linear fits.
- Omnibus RSA includes odor identity, context identity, pleasantness distance, and intensity distance. Categorical predictors remain similarity-coded (`1 = same`), so negative categorical betas indicate smaller neural distances for same-identity pairs after adjustment.
- The primary displacement coefficient model includes both rating changes plus context-pair and odor fixed effects. The context-pair-only adjustment model is retained as sensitivity.
- These analyses use 5,000 global neural condition-correspondence permutations per subject, applying one random bijection to both RDM axes, shared across ROIs and analyses and independent across subjects. Predictors remain fixed and models are refitted. This explicit RDM-level null assumes exchangeable neural condition identities under global no association; it is not a coefficient-specific conditional null or the within-run trial-label permutation used for decoding and split-half similarity.
- Report two-sided inclusive empirical p-values with add-one correction and BH-FDR q-values across ROI × scientific predictors, separately by subject/group and by model (20 omnibus tests; 10 tests per displacement model). Scientific coefficients receive permutation inference; nuisance coefficients are saved descriptively.
- Displacement scatterplots use separate univariate REML mixed models per ROI, rating, and focal semantic context (`PERSON`, `FOOD`, `LOCATION`). Pool the three pairs involving the focal context, including `CONTROL`, for 300 observations per model. Estimate one common rating slope with pair-specific fixed intercepts, an uncorrelated subject random intercept/slope, and odor and subject–odor random intercepts. Slope tests and pointwise confidence intervals use Satterthwaite degrees of freedom; FDR covers all 30 slopes. These overlapping context subsets test rating-displacement association within each pool, not a contrast against excluded-pair slopes. Display population lines equally averaging the three pair intercepts. These models have no permutation group test and do not adjust for the other rating. Preserve and report convergence and covariance-boundary diagnostics. Use `scripts/run_rating_displacement_mixed_models.m` to update only these models from saved results.
- Outputs are in `RDMs/omnibus_RSA/` and `RDMs/rating_displacement_RSA/`, each with subject results, tables, figures, combined MAT results, and a methods README. The original exploratory outputs remain in `RDMs/exploratory_RSA_results/`.

### Restricted olfactory-network searchlight

- Participant-level searchlights are complete for subjects 2-6 using the seven-mask extended support set, sniff-aligned physiology-regressed betas, semantic trials only, 4-mm physical spheres, at least 10 restricted voxels, 200 session-balanced splits, and 5,000 within-run permutations.
- `PERSON`, `FOOD`, and `LOCATION` maps are corrected jointly across all valid centers and all three contexts with a participant-level studentized maximum statistic.
- Searchlight final restriction masks, valid-center masks, effect maps, uncorrected p maps, max-stat FWE p maps, and full result/checkpoint files are in `MRI/subj_N/nifti/olf_context_similarity_searchlight/`.
- This analysis currently has participant-level inference only; it is not a group searchlight result until a prespecified common-space group model is added.
- `scripts/olf_searchlight_context_similarity_plot.m` creates an exploratory flattened ROI-topology display. It is deliberately non-anatomical and must not be presented as a brain-coordinate map.

### Global context-similarity searchlight and native T1 projections

- Completed global results for subjects 2–6 are in `MRI/subj_N/nifti/olf_context_similarity_searchlight/global/`. Saved settings confirm 6-mm spheres, at least 10 usable features, 200 session-balanced splits, and 5,000 within-run permutations, with participant-level max-stat FWE correction across centers × three semantic contexts.
- Global support is the acquired functional grid intersected with GM, positive Odor > Rest uncorrected `p < .001`, and finite GLMsingle rows; it does not imply unrestricted whole-brain coverage.
- `scripts/project_global_searchlight_to_T1.py` projects the saved maps using each participant's established `coreg/func2T1mat` and `anat/betT1brain.nii.gz` reference. Completed projections are in each global folder's `native_T1/` subdirectory. Native functional outputs are preserved.
- All maps use nearest-neighbour resampling to retain source values and significance decisions. Each subject has 19 projected maps plus six FWE-masked delta-z/null-z overlays. Background is zero for effect maps and one for p maps; the transformed valid-center mask distinguishes unsupported voxels from measured zeros. These are display-space derivatives, not new statistical tests.
- Each projection folder contains `README.md` and `projection_manifest.json` with source/reference/transform hashes, interpolation, and validation records. The runner refuses to overwrite an existing projection folder.

## Where to find group outputs

```text
MRI/group/roi_decoding_inventory/
MRI/group/context_split_half_similarity_physio/
MRI/group/roi_odor_context_template_loro/
MRI/group/roi_cross_odor_context_template_loro/
RDMs/omnibus_RSA/
RDMs/rating_displacement_RSA/
```

Large MRI and derivative files are intentionally ignored by Git. Analysis code and documentation are versioned; the existence of a local derivative should not be inferred from Git history alone.

## Checklist for future analyses

Before running or interpreting a new analysis, record:

1. Subject IDs and exact number of included trials/runs.
2. Event alignment (`sniff`, `countdown`, or `cue`) and GLMsingle derivative path.
3. Trial labels and whether `CONTROL` is included.
4. ROI selection and exact mask basenames; distinguish the focused five-ROI set from the seven-mask searchlight support set.
5. GM mask and exact functional restriction threshold (`uncorrected_p001` or `FWE_p001`).
6. Final feature counts after all intersections and finite-data checks.
7. Normalization, similarity/decoder, fold grouping, and chance/null definition.
8. Permutation exchangeability unit, number of permutations, random seed, and multiple-comparison family.
9. Subject-native versus common-space output and any transform used.
10. Output directory, production/smoke-test status, and enough saved metadata to reproduce the run.
