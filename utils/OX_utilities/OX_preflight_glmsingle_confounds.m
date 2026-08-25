function report = OX_preflight_glmsingle_confounds(varargin)
%OX_PREFLIGHT_GLMSINGLE_CONFOUNDS Audit source coverage without saving confounds.
%   REPORT = OX_PREFLIGHT_GLMSINGLE_CONFOUNDS('Subjects',2:6) checks BM
%   ordinal identity, functional/QC nframes, MRI TTL counts/regularity, and
%   raw-waveform coverage for every run. It returns all problems together.

p = inputParser;
p.addParameter('Subjects', 2:6, @(x) isnumeric(x) && isvector(x));
p.addParameter('ProjectRoot', localProjectRoot(), @(x) ischar(x) || isstring(x));
p.parse(varargin{:});
subjects = p.Results.Subjects(:)';
projectRoot = char(p.Results.ProjectRoot);
labchartDir = fullfile(projectRoot, 'labchart', 'extracted_events');
fs = 1000;
samplesPerTR = 760;
ttlThreshold = 1;

rows = cell(0, 14);
for subjidx = subjects
    subject = sprintf('subj_%d', subjidx);
    runs = OX_discover_functional_runs(subjidx, 'ProjectRoot', projectRoot, ...
        'ExpectedRuns', 80);
    bmFile = fullfile(labchartDir, sprintf('subj%d_events_bm.mat', subjidx));
    bm = load(bmFile);
    rawFile = fullfile(labchartDir, bm.bm_processing.source_raw_file);
    raw = load(rawFile, 'eventdata');
    qcFile = fullfile(projectRoot, 'motion_qc', sprintf('%s_motion_qc.mat', subject));
    q = load(qcFile, 'results');
    for runIndex = 1:numel(runs)
        nframes = bm.nframes(runIndex);
        runData = raw.eventdata{runIndex};
        rises = find(diff([false, runData(3, :) > ttlThreshold]) == 1);
        bmOnsets = bm.bm_matched_inhale_samples(:, runIndex) - ...
            bm.event_onsets(:, runIndex) .* fs;
        bmOnset = round(bmOnsets(1));
        bmIdentityOK = all(isfinite(bmOnsets)) && ...
            max(abs(bmOnsets - bmOnset)) <= 1e-6;
        eventRises = find(diff([false, runData(2, :) > 0.1]) == 1);
        bmIdentityOK = bmIdentityOK && numel(eventRises) == 23 && ...
            isequal(eventRises(5:2:end)', bm.bm_sniff_ttl_samples(:, runIndex));
        regular = numel(rises) > 1 && all(abs(diff(rises) - samplesPerTR) <= 2);
        individualOK = regular && ismember(numel(rises), [nframes, nframes + 1]);
        if individualOK && numel(rises) == nframes + 1
            frameEdges = rises(1:(nframes + 1));
        elseif individualOK && numel(rises) == nframes
            frameEdges = [rises(:); rises(end) + samplesPerTR]';
        else
            frameEdges = bmOnset + (0:nframes) .* samplesPerTR;
        end
        requiredEnd = frameEdges(end) - 1;
        coverageMargin = size(runData, 2) - requiredEnd;
        finalReplacement = coverageMargin < 0 && ...
            frameEdges(end - 1) <= size(runData, 2) + 1;
        functionalFrames = niftiFrames(runs(runIndex).functional_file);
        lengthsOK = functionalFrames == nframes && ...
            q.results.n_volumes_per_run(runIndex) == nframes && ...
            numel(q.results.fd_mm{runIndex}) == nframes && ...
            numel(q.results.robust_dvars_z{runIndex}) == nframes;
        if ~bmIdentityOK
            status = 'BM/raw identity mismatch';
        elseif ~lengthsOK
            status = 'functional/QC/nframes mismatch';
        elseif coverageMargin < 0 && ~finalReplacement
            status = 'incomplete raw respiratory coverage';
        else
            status = 'ok';
        end
        if individualOK
            alignment = 'individual MRI TTL edges';
        else
            alignment = 'BM first onset + nominal TR';
        end
        rows(end + 1, :) = {subjidx, runIndex, runs(runIndex).session, ...
            runs(runIndex).run, string(runs(runIndex).id), nframes, ...
            size(runData, 2), numel(rises), bmOnset, coverageMargin, ...
            individualOK, finalReplacement, string(alignment), string(status)}; %#ok<AGROW>
    end
    clear raw bm q
end
report = cell2table(rows, 'VariableNames', {'subject', 'run_ordinal', ...
    'session_id', 'run_id', 'run_label', 'nframes', 'raw_samples', ...
    'mri_ttl_rises', 'bm_mri_onset_sample', 'coverage_margin_samples', ...
    'individual_ttl_validated', 'final_volume_replacement_required', ...
    'alignment_method', 'status'});
end

function n = niftiFrames(filename)
info = niftiinfo(filename);
n = info.ImageSize(4);
end

function root = localProjectRoot()
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end
