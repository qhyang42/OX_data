function output = OX_rsa_context_displacement_lme(T,out,makeplots)
% Pool all three incident pairs for each semantic context, one rating slope.
T.subject=categorical(T.subject_id); T.odor=categorical(T.odor_id);
T.subject_odor=categorical(string(T.subject_id)+"_"+string(T.odor_id));
pairs=unique(string(T.context_pair),'stable'); endpoints=split(string(T.context_pair),'_vs_');
assert(size(endpoints,2)==2);
contexts=["PERSON","FOOD","LOCATION"]; rois=unique(T.roi,'stable');
predictors=["pleasantness","intensity"];
slopes=table(); fixed=table(); diagnostics=table(); residual_table=table(); predictions=table(); membership=table();
models=cell(5,2,3); convergence=cell(5,2,3); covariance=cell(5,2,3);
formula='neural_distance ~ context_pair + x + (1+x|subject) + (1|odor) + (1|subject_odor)';
for c=1:3
    include=any(endpoints==contexts(c),2);
    included_pairs=unique(string(T.context_pair(include)),'stable'); assert(numel(included_pairs)==3);
    membership=[membership;table(repmat(contexts(c),3,1),included_pairs, ...
        'VariableNames',{'focal_context','context_pair'})]; %#ok<AGROW>
    for j=1:2
        for r=1:5
            S=T(include & T.roi==rois(r),:); S.context_pair=categorical(S.context_pair,included_pairs);
            S.x=S.("abs_delta_"+predictors(j));
            assert(height(S)==300 && numel(unique(S.subject))==5 && all(isfinite(S.x)));
            assert(all(groupcounts(S.subject)==60) && all(groupcounts(S.context_pair)==100));
            assert(all(groupcounts(S.subject_odor)==3),'Each subject-odor must contribute exactly three pairs.');
            lastwarn('');
            M=fitlme(S,formula,'FitMethod','REML','CovariancePattern',{'Diagonal','Isotropic','Isotropic'},'CheckHessian',true);
            [wm,wi]=lastwarn; models{r,j,c}=M;
            [b,~,stats]=fixedEffects(M,'DFMethod','satterthwaite'); names=string(M.CoefficientNames(:));
            assert(numel(b)==4 && all(isfinite(b)) && all(isfinite(stats.pValue)));
            V=M.CoefficientCovariance; H=double(names=="x")'; assert(nnz(H)==1);
            slope=H*b; se=sqrt(H*V*H'); [pv,~,~,df]=coefTest(M,H,0,'DFMethod','satterthwaite');
            assert(all(isfinite([slope,se,pv,df])) && df>0);
            ci=slope+[-1 1]*tinv(.975,df)*se; scale=std(S.x)/std(S.neural_distance);
            [psi,mse]=covarianceParameters(M); variances=[diag(psi{1});psi{2};psi{3}];
            boundary=any(variances<max(1e-8,1e-6*mse)); status="ok";
            if boundary, status="boundary_variance"; end
            if strlength(string(wm))>0, status=status+"|fit_warning"; end
            conv=struct('fit_returned',true,'warning_id',string(wi),'warning_message',string(wm));
            if isprop(M,'ConvergenceInfo'), conv.optimizer=M.ConvergenceInfo; end
            convergence{r,j,c}=conv; covariance{r,j,c}=struct('psi',{psi},'mse',mse);
            slopes=[slopes;table(rois(r),predictors(j),contexts(c),300,5,3,slope,slope*scale,se,df, ...
                ci(1),ci(2),se*scale,ci(1)*scale,ci(2)*scale,pv,status, ...
                'VariableNames',{'roi','predictor','focal_context','n_observations','n_subjects','n_context_pairs', ...
                'slope_raw','beta_standardized','se_raw','df','ci_lower_raw','ci_upper_raw', ...
                'se_standardized','ci_lower_standardized','ci_upper_standardized','p_value','model_status'})]; %#ok<AGROW>
            fitted_y=fitted(M); residual=S.neural_distance-fitted_y;
            diagnostics=[diagnostics;table(rois(r),predictors(j),contexts(c),300,M.LogLikelihood,mse, ...
                variances(1),variances(2),variances(3),variances(4),boundary,string(wi),string(wm),status, ...
                sqrt(mean(residual.^2)),mean(residual), ...
                'VariableNames',{'roi','predictor','focal_context','n_observations','log_likelihood','residual_variance', ...
                'subject_intercept_variance','subject_slope_variance','odor_variance','subject_odor_variance', ...
                'boundary_variance','warning_id','warning_message','status','conditional_rmse','mean_residual'})]; %#ok<AGROW>
            residual_table=[residual_table;table(S.subject_id,S.odor_id,repmat(rois(r),300,1), ...
                repmat(predictors(j),300,1),repmat(contexts(c),300,1),string(S.context_pair),S.x,S.neural_distance, ...
                fitted_y,residual,residual/sqrt(mse),'VariableNames',{'subject_id','odor_id','roi','predictor', ...
                'focal_context','context_pair','absolute_rating_change','neural_distance','conditional_fitted','residual','scaled_residual'})]; %#ok<AGROW>
            fixed=[fixed;table(repmat(rois(r),numel(b),1),repmat(predictors(j),numel(b),1), ...
                repmat(contexts(c),numel(b),1),names,b,stats.SE,stats.DF,stats.pValue,stats.Lower,stats.Upper, ...
                'VariableNames',{'roi','predictor','focal_context','term','coefficient','se','df','p_value','ci_lower','ci_upper'})]; %#ok<AGROW>
            % Equal-weight marginalization of the three pair intercepts.
            xx=linspace(min(S.x),max(S.x),100)'; G=zeros(100,numel(b));
            G(:,names=="(Intercept)")=1; G(:,names=="x")=xx;
            G(:,startsWith(names,"context_pair_"))=1/3;
            yhat=G*b; yci=zeros(100,2);
            for a=1:100
                [~,~,~,df_prediction]=coefTest(M,G(a,:),0,'DFMethod','satterthwaite');
                yci(a,:)=yhat(a)+[-1 1]*tinv(.975,df_prediction)*sqrt(G(a,:)*V*G(a,:)');
            end
            assert(all(isfinite(yci(:))));
            for a=[1 100]
                grid=S([find(S.context_pair==included_pairs(1),1);find(S.context_pair==included_pairs(2),1);find(S.context_pair==included_pairs(3),1)],:);
                grid.x(:)=xx(a); yp=predict(M,grid,'Conditional',false);
                assert(abs(mean(yp)-yhat(a))<1e-10,'Marginal prediction mismatch.');
            end
            assert(abs((yhat(end)-yhat(1))/(xx(end)-xx(1))-slope)<1e-10);
            predictions=[predictions;table(repmat(rois(r),100,1),repmat(predictors(j),100,1), ...
                repmat(contexts(c),100,1),xx,yhat,yci(:,1),yci(:,2), ...
                'VariableNames',{'roi','predictor','focal_context','absolute_rating_change','prediction','ci_lower','ci_upper'})]; %#ok<AGROW>
            fprintf('[Context mixed model] %s %s %s: %s\n',contexts(c),rois(r),predictors(j),status);
        end
    end
end
assert(height(slopes)==30); slopes.q_value=OX_rsa_bh(slopes.p_value);
output=struct('grouping',"semantic_context_vs_rest",'formula',formula,'contexts',contexts,'roi_names',rois, ...
    'predictors',predictors,'slopes',slopes,'fixed_effects',fixed,'diagnostics',diagnostics, ...
    'residuals',residual_table,'predictions',predictions,'membership',membership, ...
    'models',{models},'convergence',{convergence},'covariance',{covariance});
writetable(slopes,fullfile(out,'tables','mixed_model_context_slopes.csv'));
writetable(membership,fullfile(out,'tables','mixed_model_context_membership.csv'));
writetable(fixed,fullfile(out,'tables','mixed_model_fixed_effects.csv'));
writetable(diagnostics,fullfile(out,'tables','mixed_model_diagnostics.csv'));
writetable(residual_table,fullfile(out,'tables','mixed_model_residuals.csv'));
writetable(predictions,fullfile(out,'tables','mixed_model_predictions.csv'));
fid=fopen(fullfile(out,'tables','mixed_model_context_slopes.md'),'w'); assert(fid>=0); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,['# Semantic-context versus rest: rating-displacement slopes\n\n' ...
    'Each focal context pools its three pairs, including CONTROL (300 observations, five subjects). ' ...
    'One common rating slope is fitted, with pair-specific intercepts. Separate univariate REML models for each rating and ROI; ' ...
    'two-sided Satterthwaite slope tests and BH-FDR across 30 tests. Semantic-context pools overlap and are not independent cohorts. ' ...
    'This tests rating-displacement association within each context pool, not a difference from slopes in the excluded pairs. ' ...
    'Standardization uses pooled SD(x)/SD(y) within each fitted subset. CIs are pointwise.\n\n']);
for c=contexts
    fprintf(fid,'## %s versus rest\n\n| ROI | Rating | Raw slope [95%% CI] | Standardized beta | p | q | Status |\n',c);
    fprintf(fid,'| --- | --- | --- | ---: | ---: | ---: | --- |\n');
    for i=find(slopes.focal_context==c)'
        a=slopes(i,:); fprintf(fid,'| %s | %s | %.4f [%.4f, %.4f] | %.4f | %.4g | %.4g | %s |\n', ...
            a.roi,a.predictor,a.slope_raw,a.ci_lower_raw,a.ci_upper_raw,a.beta_standardized,a.p_value,a.q_value,a.model_status);
    end
    fprintf(fid,'\n');
end
save(fullfile(out,'mixed_models.mat'),'output','-v7.3');
if makeplots, OX_plot_context_displacement_lme(output,out); end
end
