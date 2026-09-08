function coefficients = run_exploratory_rsa_analysis6_pair_slopes(varargin)
% Descriptive pooled, univariate context-pair fits for analysis 6.
% Can run independently from the saved observation CSV. Simple distance only.
p = inputParser;
p.addParameter('OutputDir', '', @(x) ischar(x) || isstring(x));
p.addParameter('MakePlots', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
root = string(p.Results.OutputDir);
if strlength(root) == 0
    root = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
        'RDMs', 'exploratory_RSA_results');
end
T = readtable(fullfile(root,'tables','analysis6_same_odor_observations.csv'), ...
    'TextType','string');
T = T(T.neural_metric == "simple",:);
assert(~isempty(T), 'No simple-distance observations found.');
rois = unique(T.roi,'stable'); pairs = unique(T.context_pair,'stable');
colors = lines(numel(pairs));
rows = struct([]);
for predictor = ["pleasantness", "intensity"]
    if p.Results.MakePlots
        f = figure('Visible','off','Color','w','Position',[50 50 1450 850]);
        layout = tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact');
    end
    for r = 1:numel(rois)
        if p.Results.MakePlots, ax = nexttile(layout); hold(ax,'on'); end
        for k = 1:numel(pairs)
            S = T(T.roi == rois(r) & T.context_pair == pairs(k),:);
            x = S.("abs_delta_" + predictor); y = S.neural_distance;
            valid = isfinite(x) & isfinite(y);
            x = x(valid); y = y(valid);
            intercept = NaN; slope = NaN; beta = NaN; status = "insufficient_data";
            if numel(x) >= 3 && std(x) > 0 && std(y) > 0
                fit = [ones(size(x)), x] \ y;
                intercept = fit(1); slope = fit(2);
                beta = slope * std(x) / std(y); status = "ok";
            end
            row = struct('roi',rois(r),'context_pair',pairs(k), ...
                'predictor',predictor,'neural_metric',"simple", ...
                'model',"pooled_univariate_linear",'n_valid',numel(x), ...
                'n_dropped',sum(~valid),'n_subjects',numel(unique(S.subject_id(valid))), ...
                'intercept',intercept,'slope_raw',slope,'beta_standardized',beta,'status',status);
            if isempty(rows), rows = row; else, rows(end+1) = row; end %#ok<AGROW>
            if p.Results.MakePlots
                scatter(ax,x,y,12,colors(k,:),'filled','MarkerFaceAlpha',.3, ...
                    'HandleVisibility','off');
                if status == "ok"
                    xx = [min(x),max(x)];
                    plot(ax,xx,intercept+slope*xx,'Color',colors(k,:),'LineWidth',2);
                end
            end
        end
        if p.Results.MakePlots
            xlabel(ax,"Absolute change in " + predictor); ylabel(ax,'Neural distance (1-r)');
            title(ax,rois(r)); yline(ax,0,':k');
        end
    end
    if p.Results.MakePlots
        ax = nexttile(layout,6); axis(ax,'off'); hold(ax,'on');
        h = gobjects(numel(pairs),1);
        for k=1:numel(pairs)
            h(k)=plot(ax,nan,nan,'Color',colors(k,:),'LineWidth',2);
        end
        legend(ax,h,replace(pairs,'_',' '),'Location','northwest');
        title(layout,["Displacement vs " + predictor + " change | simple distance"; ...
            "Separate pooled linear fits by context pair; descriptive, unadjusted"]);
        exportgraphics(f,fullfile(root,'figures','analysis6', ...
            "displacement_vs_"+predictor+"_simple.png"),'Resolution',180);
        close(f);
    end
end
coefficients = struct2table(rows);
writetable(coefficients,fullfile(root,'tables','analysis6_context_pair_slopes.csv'));
save(fullfile(root,'tables','analysis6_context_pair_slopes.mat'),'coefficients');
fid = fopen(fullfile(root,'tables','analysis6_context_pair_slopes.md'),'w');
assert(fid >= 0, 'Cannot write slope summary.');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'# Analysis 6: context-pair slopes\n\n');
fprintf(fid,['Simple neural distance (1-r), linear fits without rank transformation. ' ...
    'Each fit pools subjects within one ROI and context pair: distance = intercept + slope * absolute rating change. ' ...
    'Pleasantness and intensity are fitted separately, without covariate or subject adjustment. ' ...
    'Raw slopes match the plotted lines; standardized beta = slope * SD(x) / SD(y). ' ...
    'These are descriptive associations, not group inferential estimates.\n\n']);
fprintf(fid,'| ROI | Context pair | Pleasantness slope | Pleasantness beta | Intensity slope | Intensity beta |\n');
fprintf(fid,'| --- | --- | ---: | ---: | ---: | ---: |\n');
for r=1:numel(rois)
    for k=1:numel(pairs)
        A=coefficients(coefficients.roi==rois(r) & coefficients.context_pair==pairs(k) & coefficients.predictor=="pleasantness",:);
        B=coefficients(coefficients.roi==rois(r) & coefficients.context_pair==pairs(k) & coefficients.predictor=="intensity",:);
        fprintf(fid,'| %s | %s | %.4f | %.4f | %.4f | %.4f |\n', ...
            rois(r),pairs(k),A.slope_raw,A.beta_standardized,B.slope_raw,B.beta_standardized);
    end
end
end
