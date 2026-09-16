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
- Revised descriptive behavior analysis (2026-09-15): `scripts/behavior_group_and_displacement.m` writes to `results/behavior/group_consistency_and_trial_absolute_displacement_20260915/` (timestamp suffix on rerun). Group consistency retains each participant’s median saved Fisher-z resampling correlation with equal participant weighting. Displacement is now the absolute difference between each semantic trial rating and that participant’s same-odor CONTROL mean, taken before averaging. Individual context plots show valid trial dots; group context plots show five participant means with equal group weights. All valid trials within a context receive equal weight. Pleasantness and intensity group displacement figures share y-axis limits; group consistency figures also share limits across the two measures. Signed and odor-specific displacement analyses are removed from the active output; the superseded version is under `results/behavior/archive/`. All 4,000 trial records and missing-rating counts are retained, with 6,000 semantic trial-by-measure displacement rows including missing values. Descriptive only; PNG/PDF/FIG, CSV/MAT results, source and methods are saved.
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

- TU-only overall rating-displacement RSA is saved in `RDMs/rating_displacement_RSA_TU/`; entry point `scripts/run_rating_displacement_rsa_TU.m` reuses `RDMs/neural_TU/`. Subjects 2, 3, 4, 6 have 23, 12, 11, 39 usable features; subject 5 is explicitly insufficient (6 < 10). Both overall adjustment models use 5,000 permutations matching the TU omnibus mappings, with BH-FDR across two TU rating predictors separately per model. Neither model has q < .05. Regular context-pair results: pleasantness beta = .0422, q = .7750; intensity beta = -.0080, q = .8706. The automatic context-specific mixed-model block in `OX_formal_rsa` is commented out; the standard displacement runner now estimates overall betas only. Established context-specific results and the standalone mixed-model runner remain available.
- TU-only omnibus RSA is saved in `RDMs/omnibus_RSA_TU/`, with new TU RDMs in `RDMs/neural_TU/` and entry point `scripts/run_omnibus_rsa_TU.m`. It uses subjects 2, 3, 4, and 6; subject 5 has 6 usable voxels, below the unchanged 10-voxel minimum. The 5,000 permutations match the original mappings for included subjects. BH-FDR is across the four TU predictors; no predictor passes q < .05 (context beta = -0.0211, p = .0252, q = .1008). Existing ROI inputs/results were preserved.
- Entry points are `scripts/run_omnibus_rsa.m` and `scripts/run_rating_displacement_rsa.m`; shared implementation is in `utils/OX_utilities/`. Both reuse saved neural and behavioral RDMs for subjects 2–6 and the five olfactory ROIs, including all four contexts, with simple neural distance and linear fits.
- Omnibus RSA includes odor identity, context identity, pleasantness distance, and intensity distance. Categorical predictors remain similarity-coded (`1 = same`), so negative categorical betas indicate smaller neural distances for same-identity pairs after adjustment.
- The primary displacement coefficient model includes both rating changes plus context-pair and odor fixed effects. The context-pair-only adjustment model is retained as sensitivity.
- These analyses use 5,000 global neural condition-correspondence permutations per subject, applying one random bijection to both RDM axes, shared across ROIs and analyses and independent across subjects. Predictors remain fixed and models are refitted. This explicit RDM-level null assumes exchangeable neural condition identities under global no association; it is not a coefficient-specific conditional null or the within-run trial-label permutation used for decoding and split-half similarity.
- Report two-sided inclusive empirical p-values with add-one correction and BH-FDR q-values across ROI × scientific predictors, separately by subject/group and by model (20 omnibus tests; 10 tests per displacement model). Scientific coefficients receive permutation inference; nuisance coefficients are saved descriptively.
- Displacement scatterplots use separate univariate REML mixed models per ROI, rating, and focal semantic context (`PERSON`, `FOOD`, `LOCATION`). Pool the three pairs involving the focal context, including `CONTROL`, for 300 observations per model. Estimate one common rating slope with pair-specific fixed intercepts, an uncorrelated subject random intercept/slope, and odor and subject–odor random intercepts. Slope tests and pointwise confidence intervals use Satterthwaite degrees of freedom; FDR covers all 30 slopes. These overlapping context subsets test rating-displacement association within each pool, not a contrast against excluded-pair slopes. Display population lines equally averaging the three pair intercepts. These models have no permutation group test and do not adjust for the other rating. Preserve and report convergence and covariance-boundary diagnostics. Use `scripts/run_rating_displacement_mixed_models.m` to update only these models from saved results.
- Outputs are in `RDMs/omnibus_RSA/` and `RDMs/rating_displacement_RSA/`, each with subject results, tables, figures, combined MAT results, and a methods README. The original exploratory outputs remain in `RDMs/exploratory_RSA_results/`.
- Formal overall displacement coefficient display for the context-pair-only adjustment model is in `results/rating_displacement_RSA/` (PNG/vector PDF, subject dots, and saved group BH-FDR significance stars); reproduce with `scripts/plot_rating_displacement_RSA_formal.m`. This display reuses saved production fits and does not change the primary model.
- Current regular context-pair displacement plots are separate pleasantness and intensity figures in `results/rating_displacement_RSA_complete/separate_ratings/`; reproduce with `scripts/plot_rating_displacement_RSA_complete.m`. Group BH-FDR is recalculated across eight ROIs separately for each rating (two eight-test families), using saved permutation p-values. Coefficient models still include both ratings plus context-pair adjustment. Order: TU, AON, PirF, PirT, olfAMG, olfOFC, Putamen, A1. Controls have gray boxes/dots; TU retains n=4 and others n=5. The earlier combined plot and its original correction families remain in the parent directory for provenance.

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

### Global searchlight MNI participant-overlap displays

- `scripts/group_global_searchlight_overlap_mni.py` registers the saved global searchlight binary joint max-stat FWE `p < .05` masks for subjects 2–6 using existing `inverse(std2T1mat) × func2T1mat` transforms and nearest-neighbour interpolation to the MNI152 1-mm reference.
- `results/searchlight_context_split_half_similarity_results/` contains separate PERSON, FOOD and LOCATION count maps (0–5 significant participants), valid-center coverage counts, individual registered masks, provenance, and matching coronal PNGs plus a combined overlay. The shared montage script is `scripts/plot_odor_overlap_committee.py --analysis searchlight`.
- These are descriptive overlaps of participant-level significance, not common-space group inference. The global searchlight's acquisition and functional restriction remain applicable; display interpolation introduces no new tests.

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

### Anatomical control ROI additions (2026-09-12)

Two bilateral anatomical controls are saved in `ROIs/control_rois/`: `control_GM_putamen` (Harvard–Oxford putamen GM) and `control_A1` (Juelich TE1.0/1.1/1.2 intersected with atlas cortical GM). Each is 814 mm³ in MNI152 1-mm space, matching the rounded median 813.5 mm³ of bilateral TU, AON, PirF, PirT, olfAMG, and olfOFC; each hemisphere contributes 407 voxels. Compact patches were selected anatomically; A1 seed positioning also uses shared inverse-mapped native GM coverage, without task-response selection. A1 is feasible after this coverage positioning, so V1 was not used.

Native versions use the established affine chain and trilinear/0.2 convention and are in `MRI/subj_N/nifti/coreg/roi_decoding/control/` for subjects 2–6. This is an additional directory, not a change to primary/secondary sets or existing selection defaults. Standard-space volume matching does not imply equal native volumes or usable features. Both controls have ≥10 finite sniff GLMsingle features after native GM intersection in all five participants; several Odor > Rest restricted intersections are insufficient. No existing functional restriction default is changed. See `ROIs/control_rois/README.md`, `manifest.json`, `native_voxel_counts.csv`, `feature_qc.csv`, and the QC overlays for definitions, counts, provenance and limitations. Producing scripts: `scripts/build_control_rois.py` and `scripts/check_control_rois.py`.


White-matter addition and control-specific restriction update (2026-09-12): `control_WM` adds an 814-mm³ bilateral deep cerebral WM control (407 voxels/hemisphere), from Harvard–Oxford WM probability ≥90% with a 3-mm interior margin. MNI/native files share the above control directories; `white_matter_manifest.json`, `white_matter_native_qc.csv`, and `white_matter_QC.png` record provenance and QA. Native anatomical counts for subjects 2–6 are 137/134/147/24/108. Per explicit user instruction, **all control ROIs omit native GM intersection and functional restriction**; prior restriction-QC statuses do not determine control eligibility. Existing non-control analysis defaults are unchanged. Saved putamen/A1 masks are unchanged. Canonical GM-only GLMsingle outputs cannot supply full WM betas; WM analysis requires inputs including WM voxels. The mask-creation task does not change generic analysis runners.

### Putamen/A1 control context split-half analysis (2026-09-12)

Production results are in `MRI/group/context_split_half_similarity_physio_control/`. Subjects 2–6, putamen and A1 only; canonical sniff-aligned physio GLMsingle inputs, **native GM intersection but no functional restriction**, as explicitly requested for this analysis. This GM-only restriction is an analysis-specific update to the earlier control policy. Uses the established semantic-only within-run centering, 600 semantic trials, 80 exact runs, 20 odors, 200 session-balanced splits, 5,000 complete within-run label permutations and equal-weight subject-null aggregation. The original olfactory result folder is unchanged.

Final usable putamen/A1 features by subject 2–6: 59/78, 37/95, 46/66, 107/42, 103/13. All ten subject–ROI analyses completed. Overall mean delta-z: putamen −0.06031 (p=0.84083, BH q=0.84083); A1 +0.10819 (p=0.22416, BH q=0.44831). Neither overall control effect is significant. No ROI-by-context group effect survives the separate six-test BH-FDR family. Overall correction is across two control ROIs; context-specific correction across six control ROI × context tests, separately from olfactory results. Inference pertains to the measured participants, not population random effects. Participant effects and nulls remain in the MAT/CSV outputs.

Runner: `scripts/context_split_half_similarity_control_analysis.m`; FDR/provenance: `scripts/summarize_control_split_half.py`. `OX_roi_context_split_half_similarity` now accepts `UseFunctionalRestriction` (default true); the control runner explicitly sets false and uses the existing flat-directory override to select `roi_decoding/control/`. No legacy masks or white-matter ROI enter this analysis.

### Putamen/A1 control omnibus RSA (2026-09-12)

Completed in `RDMs/omnibus_RSA_control/`, with new neural RDMs in `RDMs/neural_control/`; runner `scripts/run_omnibus_rsa_control.m`. Includes subjects 2–6 and native putamen/A1 masks intersected with GM and finite betas, **without functional restriction**. Final usable counts match the control split-half analysis. Unlike the semantic-only split-half model, omnibus uses all 800 trials and all four contexts for within-run centering and 80-condition RDM construction. Pearson distance only, linear standardized omnibus model with odor/context identity (1=same), pleasantness and intensity distances.

All 5,000 global neural condition-correspondence permutations match the original omnibus mappings. Equal-weight five-subject null aggregation and two-sided inclusive add-one p-values; BH-FDR across **eight control ROI × scientific predictor tests**, separately per subject and at group level, separate from prior olfactory results. No group coefficient passes q < .05. A1 context beta = −0.022583 (p=.017596, q=.115177); A1 intensity beta = +0.013255 (p=.028794, q=.115177). These are not corrected-significant effects or population random-effects estimates. Signed subject coefficients/nulls and absolute-beta display figures are retained.

The shared RDM builder now accepts optional `UseFunctionalRestriction=false` and `ComputeCrossnobis=false`, with both defaults remaining true. This control run disables crossnobis because it is not used in omnibus RSA; an initial auxiliary whitening failure for subj4 A1 is recorded separately in `tmp/omnibus_control_crossnobis_attempt/`, with no ROI/subject exclusion or threshold relaxation. New simple RDMs match the earlier full-builder simple RDMs exactly for subjects 2/3 and match independent Pearson calculations for all subjects. Existing olfactory and TU outputs are preserved.

### Putamen/A1 control rating-displacement RSA (2026-09-14)

Completed in `RDMs/rating_displacement_RSA_control/`, using saved simple-distance `RDMs/neural_control/` inputs (native GM, no functional restriction), all five subjects, all four contexts and 120 same-odor cross-context pairs per subject/ROI. Runner: `scripts/run_rating_displacement_rsa_control.m`. Primary context-pair-plus-odor adjusted and context-pair-only sensitivity models each have 5,000 global condition-correspondence permutations matching control omnibus RSA, equal-weight subject-null aggregation, two-sided inclusive add-one p-values and separate four-test control ROI × rating BH-FDR families per model/subject or group. Neither model has an FDR-significant group coefficient. Primary pleasantness/intensity betas: putamen −0.04111/−0.03205, A1 +0.06813/−0.01747 (q=.68946–.73345). Context-pair-only q=.83303 for all four group tests.

The twelve supplementary context-pooled REML models were attempted separately. Ten return finite inference values, all with boundary/covariance warnings; two putamen pleasantness fits (FOOD and LOCATION pools) fail inference and are explicitly unavailable. No simplification or exclusion was performed. `mixed_model_failures.csv`, diagnostic tables and labeled plots preserve this limitation. Mixed slope FDR retains the full twelve-test family with unavailable slots conservatively set to p=1 for adjustment only, while failed p/q remain unreported. These supplementary unstable fits do not invalidate the completed permutation RSA, and should not support firm inferential claims. Validation/provenance are saved in the result folder; prior olfactory outputs are unchanged.

### RSA-region temporal SNR control figure (2026-09-16)

Saved in `results/ROI_tSNR/` as PNG, vector PDF/SVG, subject/run CSVs, exact input metadata and voxelwise run caches. Order: TU, AON, PirF, PirT, olfAMG, olfOFC, Putamen, A1. Uses the saved RSA feature indices: olfactory GM + positive Odor > Rest uncorrected p < .001; controls GM only; finite-beta selection retained. Temporal SNR is voxelwise temporal mean/sample SD in each smoothed, realigned `sr*.nii` GLMsingle-input run, before GLMsingle nuisance regression, averaged over voxels and then equally over 80 runs. No additional temporal filtering, detrending or frame censoring. Boxplots summarize subject values with subject overlays; descriptive, no inference.

All subjects 2–6 retained except established insufficient TU subject 5 (6 < 10 RSA features). Three subject-4 A1 voxels have zero temporal SD in some runs; these are explicitly logged and removed across all runs for this tSNR summary, leaving 63 of 66 A1 RSA features. All other feature counts unchanged. This does not alter RSA inputs or results. Final tSNR feature counts are in `subject_tsnr.csv`, affected voxel/run records in `invalid_voxel_runs.csv`. Verified 400 functional runs against canonical metadata/events, 39 usable subject–ROI summaries and 3,120 run-level measurements. Reproduce with `scripts/export_rsa_tsnr_inputs.m` and `scripts/plot_rsa_roi_tsnr.py`.
