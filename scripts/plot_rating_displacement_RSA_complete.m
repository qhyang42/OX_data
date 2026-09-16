function plot_rating_displacement_RSA_complete(separate_ratings)
if nargin<1, separate_ratings=true; end
% Combine saved overall context-pair coefficients, retaining inference families.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'utils','OX_utilities'));
sources=["rating_displacement_RSA_TU","rating_displacement_RSA","rating_displacement_RSA_control"];
S=cell(1,3); M=cell(1,3);
for k=1:3
    L=load(fullfile(root,'RDMs',sources(k),'results.mat'),'results'); S{k}=L.results;
    keep=arrayfun(@(m) string(m.name)=="context_pair",S{k}.models);
    assert(nnz(keep)==1); M{k}=S{k}.models(keep);
    assert(isequal(M{k}.predictors,["pleasantness","intensity"]));
    assert(all(isfinite(M{k}.beta),'all') && all(M{k}.group_statistics.n_permutations==5000));
end
assert(isequal(S{1}.subject_ids,[2 3 4 6]) && isequal(S{2}.subject_ids,2:6) && isequal(S{3}.subject_ids,2:6));
assert(isequal(S{1}.roi_names,"olf_TU"));
assert(isequal(S{2}.roi_names,["AON","PirF","PirT","olfAMG","olfOFC"]));
assert(isequal(S{3}.roi_names,["control_GM_putamen","control_A1"]));
out=fullfile(root,'results','rating_displacement_RSA_complete');
if separate_ratings, out=fullfile(out,'separate_ratings'); end
assert(~isfolder(out),'Output exists; preserve previous figure.'); mkdir(out);
P=struct('analysis',"displacement",'subject_ids',2:6, ...
    'roi_names',[S{1}.roi_names,S{2}.roi_names,S{3}.roi_names], ...
    'roi_display_names',["TU","AON","PirF","PirT","olfAMG","olfOFC","Putamen","A1"], ...
    'control_roi_names',S{3}.roi_names,'output_dir',out,'figure_output_dir',out);
G=[M{1}.group_statistics;M{2}.group_statistics;M{3}.group_statistics];
P.models=struct('name',"context_pair",'predictors',M{1}.predictors,'beta',nan(5,8,2),'group_statistics',G);
[found,rows]=ismember(S{1}.subject_ids,P.subject_ids); assert(all(found));
P.models.beta(rows,1,:)=M{1}.beta;
P.models.beta(:,2:6,:)=M{2}.beta; P.models.beta(:,7:8,:)=M{3}.beta;
assert(all(isnan(P.models.beta(4,1,:)),'all') && nnz(isfinite(P.models.beta))==78);
for r=1:8
    for j=1:2
        idx=G.roi==P.roi_names(r) & G.predictor==P.models.predictors(j);
        b=P.models.beta(:,r,j);
        assert(nnz(idx)==1 && nnz(isfinite(b))==G.n_subjects(idx));
        assert(abs(mean(b,'omitnan')-G.mean_beta(idx))<1e-12);
    end
end
if separate_ratings
    G.original_q_value=G.q_value;
    for predictor=["pleasantness","intensity"]
        idx=G.predictor==predictor; assert(nnz(idx)==8);
        G.q_value(idx)=OX_rsa_bh(G.p_value(idx));
    end
    P.models.group_statistics=G;
    save(fullfile(out,'plot_data.mat'),'P','sources');
    writetable(G,fullfile(out,'group_statistics.csv'));
    % Subject p-values retained; no subject-level significance is plotted.
    ST=[M{1}.subject_statistics;M{2}.subject_statistics;M{3}.subject_statistics];
    ST.Properties.VariableNames{'q_value'}='original_q_value';
    writetable(ST,fullfile(out,'subject_statistics.csv'));
    plot_separate(P);
    fid=fopen(fullfile(out,'README.md'),'w'); assert(fid>=0);
    fprintf(fid,['# Rating displacement RSA: separate rating families\n\n' ...
        'Pleasantness and intensity are displayed separately. Group BH-FDR is recomputed across all eight ROIs within each rating (two separate eight-test families), from saved two-sided permutation p-values. ' ...
        'Original q-values are retained in a separate column. Stars: * q < .05, ** q < .01, *** q < .001. ' ...
        'Coefficient estimation is unchanged: regular context-pair adjustment with both rating changes in the same model, no odor fixed effects. Separate analysis here means separate display and multiplicity families, not new univariate fits.\n\n' ...
        'Order: TU, AON, PirF, PirT, olfAMG, olfOFC, Putamen, A1. Controls use gray boxes/dots (GM only); other ROIs use GM plus positive Odor > Rest p < .001. ' ...
        'Boxes: median, IQR, whiskers to extreme values within 1.5 IQR. Dots: subjects 2–6 with fixed offsets. TU excludes subject 5 (6 usable features < 10), so n=4; all other ROIs n=5. ' ...
        'Plots share y limits. All six context pairs and 20 odors contribute. Saved source results: RDMs/rating_displacement_RSA_TU, rating_displacement_RSA, rating_displacement_RSA_control. ' ...
        '5000 global condition-correspondence permutations; equal-weight subject null aggregation; inference concerns measured subjects. No fits or permutations rerun. ' ...
        'Subject statistics retain original q-values only; new correction applies to plotted group tests.\n\n' ...
        'Reproduce: scripts/plot_rating_displacement_RSA_complete.m (default separate_ratings=true). Earlier combined figure is preserved in parent directory.\n']);
    fclose(fid);
    return
end
P.figure_note=["Gray controls: GM only | Left box: pleasantness; right box: intensity"; ...
    "BH-FDR families: original ROIs 10; TU 2; controls 4 | * q < .05, ** q < .01, *** q < .001"];
OX_plot_formal_rsa(P);
save(fullfile(out,'plot_data.mat'),'P','sources');
writetable(G,fullfile(out,'group_statistics.csv'));
writetable([M{1}.subject_statistics;M{2}.subject_statistics;M{3}.subject_statistics],fullfile(out,'subject_statistics.csv'));
fid=fopen(fullfile(out,'README.md'),'w'); assert(fid>=0);
fprintf(fid,['# Complete rating displacement RSA\n\n' ...
    'Overall signed standardized coefficients from the regular context_pair model: both rating changes plus context-pair fixed effects, no odor fixed effects. ' ...
    'Order: TU, AON, PirF, PirT, olfAMG, olfOFC, Putamen, A1. All six cross-context pairs and 20 odors contribute. ' ...
    'Saved production results are reused from RDMs/rating_displacement_RSA_TU, rating_displacement_RSA, and rating_displacement_RSA_control. No fits were rerun.\n\n' ...
    'Boxes show median and IQR, whiskers extend to extreme values within 1.5 IQR; dots show individual subjects. ' ...
    'Subjects 2–6 have fixed left-to-right offsets. TU has subjects 2, 3, 4, 6; subject 5 is explicitly missing (6 usable features, below 10). ' ...
    'All other ROIs have five subjects. Native maps are not combined; only subject coefficients are aggregated. ' ...
    'Controls have gray boxes and dots; within each ROI pleasantness is left and intensity right. ' ...
    'Controls use GM only, whereas six olfactory ROIs use GM plus positive Odor > Rest uncorrected p < .001.\n\n' ...
    'Stars preserve saved signed-beta permutation inference and separate BH-FDR families: 10 original ROI × rating tests, 2 TU tests, 4 control tests. ' ...
    'No joint 16-test correction was introduced. * q < .05, ** q < .01, *** q < .001; unmarked means q >= .05. ' ...
    'Each group null averages selected subject nulls equally across 5000 global condition-correspondence permutations. ' ...
    'Inference concerns measured subjects, not population random effects. See source READMEs for full methods and feature counts.\n\n' ...
    'Reproduce with scripts/plot_rating_displacement_RSA_complete.m. Files: PNG, vector PDF, plot_data.mat, subject and group CSVs.\n']);
fclose(fid);
end

function plot_separate(P)
M=P.models; G=M.group_statistics; colors=lines(5); pc=lines(4);
lo=min(M.beta,[],'all','omitnan'); hi=max(M.beta,[],'all','omitnan'); span=hi-lo;
for j=1:2
    pred=M.predictors(j);
    f=figure('Visible','off','Color','w','Position',[50 50 1350 650]);
    ax=axes(f); hold(ax,'on');
    for r=1:8
        b=M.beta(:,r,j); valid=isfinite(b); offsets=linspace(-.12,.12,5)';
        bc=pc(j+2,:); dc=colors(valid,:);
        if ismember(P.roi_names(r),P.control_roi_names)
            bc=[.55 .55 .55]; dc=repmat([.4 .4 .4],nnz(valid),1);
        end
        boxchart(ax,repmat(r,nnz(valid),1),b(valid),'BoxFaceColor',bc,'BoxWidth',.48,'MarkerStyle','none','HandleVisibility','off');
        scatter(ax,r+offsets(valid),b(valid),42,dc,'filled','HandleVisibility','off');
        idx=G.roi==P.roi_names(r) & G.predictor==pred; q=G.q_value(idx);
        star=""; if q<.001, star="***"; elseif q<.01, star="**"; elseif q<.05, star="*"; end
        if strlength(star)>0
            text(ax,r,max(b,[],'omitnan')+.035*span,star,'HorizontalAlignment','center','FontSize',16,'FontWeight','bold');
        end
    end
    yline(ax,0,':k','HandleVisibility','off');
    xlim(ax,[.5 8.5]); ylim(ax,[lo-.10*span hi+.18*span]);
    ax.XTick=1:8; ax.XTickLabel=P.roi_display_names+" (n="+string(sum(isfinite(M.beta(:,:,j)),1))+")";
    set(ax,'FontName','Arial','FontSize',12,'LineWidth',1,'TickDir','out','Box','off','TickLabelInterpreter','none');
    ylabel(ax,'Standardized displacement coefficient (β)');
    title(ax,["Rating displacement RSA | "+upper(extractBefore(pred,2))+extractAfter(pred,1); ...
        "Context-pair adjustment | BH-FDR across 8 ROIs within this rating"; ...
        "Gray: controls | * q < .05, ** q < .01, *** q < .001"],'FontWeight','normal','Interpreter','none');
    h=gobjects(5,1);
    for s=1:5, h(s)=scatter(ax,nan,nan,42,colors(s,:),'filled'); end
    legend(ax,h,"Subject "+string(P.subject_ids),'Location','eastoutside');
    stem=fullfile(P.output_dir,"coefficients_"+pred);
    exportgraphics(f,stem+".png",'Resolution',180);
    exportgraphics(f,stem+".pdf",'ContentType','vector'); close(f);
end
end
