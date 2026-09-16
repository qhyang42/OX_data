# Behavior: group consistency and CONTROL displacement

## Outputs

- `consistency/group/`: one group consistency figure per rating measure.
- `displacement/individual/`: each participant's absolute trial-displacement box charts, separately for intensity and pleasantness.
- `displacement/group/`: group absolute trial-displacement summaries.
- `tables/`: all trial identifiers/ratings, condition means and valid/missing counts, participant and group consistency summaries, trial-level absolute displacement and participant/group displacement summaries.
- Each figure is available as PNG, vector PDF, and editable MATLAB FIG. `behavior_results.mat` stores arrays and tables; `behavior_group_and_displacement_source.m` is the producing source snapshot.

## Inputs and checks

Subjects 2–6, 800 trials, 80 runs, 20 odors and 200 trials per context per participant. Exact labels and session/run/trial identifiers are loaded using `OX_load_trial_metadata` and verified against `behavior/subj_N/behavior.mat`. Saved consistency trial inputs are required to match current ratings and metadata exactly, including NaNs. Odor labels are verified identical across participants. No subject, run, odor or condition is excluded. All 4,000 trials are retained in the table. Missing values are omitted only from measure-specific means, plotted trial dots and the original resampling; no imputation. Condition counts can differ (subject 4 has 8–12 trials per odor/context).

Ratings are on the original preprocessed scale: screen x position minus 756 (rating units, not a standardized or 0–100 scale). Pleasantness uses `valence_all`; intensity uses `intensity_all`.

## Group consistency

Reuse the exact saved `results/behavior/subj_N_MEASURE_stats.mat` resamples from `scripts/behavior_analysis.m` (twister seed 20260909). Each context has 100 profiles, each drawing one available rating per odor; Pearson correlations across 20 odors yield 9,900 directed off-diagonal correlations. Cross-context comparisons use 9,900 sampled pairs. The mixed-context reference samples a context separately per odor to generate 100 profiles and 9,900 correlations. This is a descriptive reference, not a formal permutation null.

Each dot is one participant's median Fisher-z correlation, matching the median-centered interpretation of the original boxplots. Each black diamond is the arithmetic mean of the five participant medians, on the Fisher-z scale. The dashed line is the mean participant mixed-context median. Boxes summarize the five participant medians, not pooled correlated resamples. Participant mean z values are also saved, but do not define the plotted dots. Subject ID labels identify dots. No new random resampling is performed. Both group consistency figures use identical y-axis limits, rounded outward to tenths from the minimum and maximum participant median across both measures.

## Trial-level absolute CONTROL displacement (revised)

For each semantic-context trial, absolute displacement is `abs(trial rating - mean CONTROL rating for the same subject, odor and measure)`. The CONTROL reference uses all available CONTROL repeats for that odor; it is not a temporal pre-exposure baseline. Odor identity is used only to match the baseline and preserve provenance. No odor-level displacement summaries or plots and no signed-displacement analyses are produced.

Individual box charts show up to 200 valid trial deviations in each semantic context, with every valid trial overlaid as a dot. Missing rating trials remain as NaN in the trial-level output and are counted explicitly; each baseline must have at least one valid rating. Black diamonds show the mean absolute deviation over available trials in that context. All valid trials receive equal weight, so odors with different observed counts may contribute different weights.

Group boxes summarize five participant context means. Each labeled dot is one participant's mean absolute trial deviation; the black diamond averages those five means equally. Pleasantness and intensity group plots use identical y-axis limits starting at zero, with the upper limit calculated jointly from all participant means in both measures. Trial counts are not treated as independent group sample sizes. Group CSV SDs describe variation across participants, not SEMs.

Absolute values are taken at the trial level before any averaging. This differs from the superseded absolute difference between condition means: within-context trial variability contributes to the new measure. The same CONTROL mean is reused for all matching semantic trials. No CONTROL self-deviation box is included; CONTROL defines the reference, and the analyzed contexts are FOOD, PERSON and LOCATION.

All boxes show median and quartiles with 1.5-IQR whiskers; dots include all valid observations, including values outside whiskers. Dot offsets are deterministic for legibility.

## Interpretation

These are descriptive reanalyses. No hypothesis tests, confidence intervals, correction families, or population random-effects claims are introduced. Reused trials, odor-level observations, and resampled correlations are not treated as independent substitutes for subjects. A nonzero absolute displacement can arise from trial variability and measurement noise and is not itself evidence of a context effect.

Reproduce from the project root with MATLAB R2025b: `run('scripts/behavior_group_and_displacement.m')`. Existing output folders are preserved through a timestamp suffix. The superseded condition-mean/signed/odor-plot analysis is archived under `results/behavior/archive/`; it is not the current analysis.
