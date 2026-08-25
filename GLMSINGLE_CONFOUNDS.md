# OX GLMsingle run-wise nuisance regressors

`OX_build_glmsingle_confounds` builds a standalone derivative; it does not
modify preprocessing, functional images, registration, masks, event labels,
or existing GLMsingle outputs.

## Authoritative inputs and ordering

- `labchart/extracted_events/subjN_events_bm.mat` is the source of truth for run ordinal,
  `nframes`, MRI onset, saved sniff TTLs, and BreathMetrics provenance.
- The 1-kHz waveform is read only from the
  `labchart/extracted_events/subjN_events_raw.mat` named by
  `bm_processing.source_raw_file`. Raw cell
  ordinal is not independently reinterpreted or reordered. Its event TTL
  positions and MRI onset must exactly match the saved BM metadata.
- `OX_discover_functional_runs` supplies numeric session/run order and pairs
  each `sr*.nii` with its `rp_*.txt`.
- Saved motion QC run IDs and `nframes` must exactly match discovery and BM
  metadata.

## Respiratory processing

For every run, the builder applies the shared OX polarity map (subject 3,
runs 1-10 are multiplied by -1), rejects any nonfinite sample, and applies
MATLAB `lowpass(airflow,10,1000)`. It extracts the complete acquisition
segment and applies `detrend(x,'linear')` to remove its least-squares
intercept and linear trend. Airflow is averaged in each MRI window.

Squared airflow is computed from the cleaned 1-kHz signal before frame
averaging. Respiratory volume is computed continuously as
`cumtrapz(cleaned_airflow)/1000`, never reset at breaths, and the integrated
trace is then linearly detrended before frame averaging. The three
frame-matched columns are saved before standardization and independently
z-scored within run. No HRF convolution is applied.

Individual MRI TTL boundaries are used when the recording contains exactly
`nframes` or `nframes+1` regular rising edges (759-761 samples apart). For
`nframes` edges, only the final boundary is extrapolated by the known
760-sample TR. Otherwise the BM-authoritative first MRI onset plus the known
TR defines all windows. The waveform must cover every complete window.
The sole exception is a shortage confined to the final volume: its three
pre-z-scored respiratory values are copied from volume `N-1`, and volume `N`
is forced into the bad-volume set with an explicit `incomplete respiration`
reason. A shortage affecting any earlier volume remains fatal.

## Final matrix

Each run is:

```text
[6 motion, 6 derivatives, 6 motion squares, 6 derivative squares,
 airflow_z, airflow_sq_z, resp_volume_z,
 one unique one-hot column per (FD > 0.4 OR robust DVARS z > 5 OR
 final-volume respiratory replacement) volume]
```

This is `27 + number_of_bad_volumes` columns. Spike reasons distinguish
`FD only`, `DVARS only`, `both`, and respiratory-replacement combinations.
Volumes are never removed.

## Commands

From MATLAB after `setup_ox`:

```matlab
preflight = OX_preflight_glmsingle_confounds('Subjects', 2:6);
runConfounds = OX_build_glmsingle_confounds(3);  % one subject
run('scripts/run_glmsingle_confounds.m');        % subjects 2-6, all-or-none
```

The all-subject script stops before building if preflight detects any
incomplete waveform coverage. Once all subjects exist, compact QC is:

```matlab
OX_summarize_glmsingle_confounds('Subjects', 2:6);
```

## GLMsingle compatibility

The installed `GLMestimatesingletrial.m` documentation explicitly states
that `opt.extraregressors` may contain a different number of extra
regressors in each run. Its implementation validates the number of cells,
then passes each run's matrix independently into `GLMestimatemodel`, whose
nuisance projection is also constructed per run. No equal-column padding is
needed or appropriate.
