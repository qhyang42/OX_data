function output = OX_rsa_displacement_lme(T,out,makeplots,previous)
% Marginal context-pair slopes from univariate mixed models, no permutations.
% Legacy pair-specific analysis. The current formal pipeline uses
% OX_rsa_context_displacement_lme for pooled semantic-context versus rest fits.
% Optional previous output regenerates tables/plots without refitting models.
T.subject=categorical(T.subject_id); T.odor=categorical(T.odor_id);
T.subject_odor=categorical(string(T.subject_id)+"_"+string(T.odor_id));
pairs=unique(string(T.context_pair),'stable');
T.context_pair=categorical(T.context_pair,pairs);
rois=unique(T.roi,'stable'); colors=lines(6);
slopes=table(); fixed=table(); diagnostics=table(); residual_table=table(); predictions=table();
fitted_models=cell(5,2); convergence_details=cell(5,2); covariance_details=cell(5,2);
formula='neural_distance ~ context_pair*x + (1+x|subject) + (1|odor) + (1|subject_odor)';
for j=1:2
    predictors=["pleasantness","intensity"]; predictor=predictors(j);
    if makeplots
        f=figure('Visible','off','Color','w','Position',[50 50 1450 900]);
        layout=tiledlayout(f,2,3,'TileSpacing','compact');
        diagnostic_figure=figure('Visible','off','Color','w','Position',[50 50 1150 1500]);
        diagnostic_layout=tiledlayout(diagnostic_figure,5,2,'TileSpacing','compact');
    end
    for r=1:5
        S=T(T.roi==rois(r),:); S.x=S.("abs_delta_"+predictor);
        assert(height(S)==600 && numel(unique(S.subject))==5 && all(isfinite(S.x)));
        if nargin>=4
            M=previous.models{r,j};
            warning_message=previous.convergence{r,j}.warning_message;
            warning_id=previous.convergence{r,j}.warning_id;
        else
            lastwarn('');
            M=fitlme(S,formula,'FitMethod','REML','CovariancePattern',{'Diagonal','Isotropic','Isotropic'},'CheckHessian',true);
            [warning_message,warning_id]=lastwarn;
        end
        fitted_models{r,j}=M;
        [b,~,stats]=fixedEffects(M,'DFMethod','satterthwaite'); names=string(M.CoefficientNames(:));
        assert(all(isfinite(b)) && all(isfinite(stats.pValue)));
        V=M.CoefficientCovariance;
        [psi,mse]=covarianceParameters(M); variances=vertcat(psi{2},psi{3}); variances=[diag(psi{1});variances];
        boundary=any(variances<max(1e-8,1e-6*mse));
        status="ok"; if boundary, status="boundary_variance"; end
        if strlength(string(warning_message))>0, status=status+"|fit_warning"; end
        conv=struct('fit_returned',true,'warning_id',string(warning_id),'warning_message',string(warning_message));
        if isprop(M,'ConvergenceInfo'), conv.optimizer=M.ConvergenceInfo; end
        convergence_details{r,j}=conv; covariance_details{r,j}=struct('psi',{psi},'mse',mse);
        conditional_fitted=fitted(M); residual=S.neural_distance-conditional_fitted;
        diagnostics=[diagnostics;table(rois(r),predictor,600,M.LogLikelihood,mse, ...
            variances(1),variances(2),variances(3),variances(4),boundary,string(warning_id),string(warning_message),status, ...
            sqrt(mean(residual.^2)),mean(residual), ...
            'VariableNames',{'roi','predictor','n_observations','log_likelihood','residual_variance', ...
            'subject_intercept_variance','subject_slope_variance','odor_variance','subject_odor_variance', ...
            'boundary_variance','warning_id','warning_message','status','conditional_rmse','mean_residual'})]; %#ok<AGROW>
        residual_table=[residual_table;table(S.subject_id,S.odor_id,repmat(rois(r),600,1), ...
            repmat(predictor,600,1),string(S.context_pair),conditional_fitted,residual,residual/sqrt(mse), ...
            'VariableNames',{'subject_id','odor_id','roi','predictor','context_pair','conditional_fitted','residual','scaled_residual'})]; %#ok<AGROW>
        fixed=[fixed;table(repmat(rois(r),numel(b),1),repmat(predictor,numel(b),1),names(:), ...
            b,stats.SE,stats.DF,stats.pValue,stats.Lower,stats.Upper, ...
            'VariableNames',{'roi','predictor','term','coefficient','se','df','p_value','ci_lower','ci_upper'})]; %#ok<AGROW>
        if makeplots, ax=nexttile(layout,r); hold(ax,'on'); end
        for k=1:6
            idx=string(S.context_pair)==pairs(k); n=nnz(idx); assert(n==100);
            % Match fitted coefficient names and verify each contrast against
            % the model's own marginal predictions at x=0 and x=1.
            grid=S(repmat(find(idx,1),2,1),:); grid.x=[0;1];
            [yp,~]=predict(M,grid,'Conditional',false);
            H=zeros(1,numel(b)); H(names=="x")=1;
            if k>1
                term="context_pair_"+pairs(k)+":x";
                reverse="x:context_pair_"+pairs(k);
                which=(names==term | names==reverse); assert(nnz(which)==1);
                H(which)=1;
            end
            slope=H*b; assert(abs(diff(yp)-slope)<1e-10,'Mixed-model slope contrast mismatch.');
            se=sqrt(H*V*H'); [pv,~,~,df]=coefTest(M,H,0,'DFMethod','satterthwaite');
            assert(isfinite(df)&&isfinite(pv)&&isfinite(se));
            ci=slope+[-1 1]*tinv(.975,df)*se;
            scale=std(S.x(idx))/std(S.neural_distance(idx));
            slopes=[slopes;table(rois(r),predictor,pairs(k),n,5,slope,slope*scale,se,df,ci(1),ci(2), ...
                se*scale,ci(1)*scale,ci(2)*scale,pv,status, ...
                'VariableNames',{'roi','predictor','context_pair','n_observations','n_subjects','slope_raw', ...
                'beta_standardized','se_raw','df','ci_lower_raw','ci_upper_raw','se_standardized', ...
                'ci_lower_standardized','ci_upper_standardized','p_value','model_status'})]; %#ok<AGROW>
            grid=S(repmat(find(idx,1),100,1),:);
            grid.x=linspace(min(S.x(idx)),max(S.x(idx)),100)';
            [yhat,yci]=predict(M,grid,'Conditional',false,'Prediction','curve','DFMethod','satterthwaite');
            assert(max(abs(yhat-(yp(1)+slope*grid.x)))<1e-10);
            predictions=[predictions;table(repmat(rois(r),100,1),repmat(predictor,100,1), ...
                repmat(pairs(k),100,1),grid.x,yhat,yci(:,1),yci(:,2), ...
                'VariableNames',{'roi','predictor','context_pair','absolute_rating_change','prediction','ci_lower','ci_upper'})]; %#ok<AGROW>
            if makeplots
                fill(ax,[grid.x;flipud(grid.x)],[yci(:,1);flipud(yci(:,2))],colors(k,:), ...
                    'FaceAlpha',.10,'EdgeColor','none','HandleVisibility','off');
                scatter(ax,S.x(idx),S.neural_distance(idx),10,colors(k,:),'filled','MarkerFaceAlpha',.26,'HandleVisibility','off');
                plot(ax,grid.x,yhat,'Color',colors(k,:),'LineWidth',2,'HandleVisibility','off');
            end
        end
        if makeplots
            xlabel(ax,"Absolute change in "+predictor); ylabel(ax,'Neural distance (1-r)');
            label=rois(r);
            if strlength(string(warning_message))>0, label=label+" (covariance warning)";
            elseif boundary, label=label+" (boundary variance)"; end
            title(ax,label); yline(ax,0,':k');
            da=nexttile(diagnostic_layout,2*r-1);
            scatter(da,conditional_fitted,residual,8,'filled'); yline(da,0,':k'); title(da,rois(r));
            xlabel(da,'Conditional fitted distance'); ylabel(da,'Residual');
            da=nexttile(diagnostic_layout,2*r); qqplot(da,residual); title(da,rois(r)+" residual Q-Q");
        end
        fprintf('[Mixed model] %s %s: %s\n',rois(r),predictor,status);
    end
    if makeplots
        ax=nexttile(layout,6); hold(ax,'on'); h=gobjects(6,1);
        for k=1:6, h(k)=plot(ax,nan,nan,'Color',colors(k,:),'LineWidth',2); end
        legend(ax,h,replace(pairs,'_',' '),'Location','northwest'); axis(ax,'off');
        title(layout,["Displacement vs "+predictor+" change | simple distance"; ...
            "Univariate mixed models: context-pair slopes and population predictions with 95% confidence bands"]);
        saveplot(f,fullfile(out,'figures',"displacement_vs_"+predictor+"_mixed_model"));
        title(diagnostic_layout,predictor+" mixed-model residual diagnostics");
        saveplot(diagnostic_figure,fullfile(out,'figures',"mixed_model_diagnostics_"+predictor));
    end
end
assert(height(slopes)==60); slopes.q_value=OX_rsa_bh(slopes.p_value);
writetable(slopes,fullfile(out,'tables','mixed_model_pair_slopes.csv'));
writetable(fixed,fullfile(out,'tables','mixed_model_fixed_effects.csv'));
writetable(diagnostics,fullfile(out,'tables','mixed_model_diagnostics.csv'));
writetable(residual_table,fullfile(out,'tables','mixed_model_residuals.csv'));
writetable(predictions,fullfile(out,'tables','mixed_model_predictions.csv'));
fid=fopen(fullfile(out,'tables','mixed_model_pair_slopes.md'),'w'); assert(fid>=0); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'# Mixed-model context-pair slopes\n\nSeparate univariate REML models per ROI and rating. Satterthwaite tests; BH-FDR over all 60 slopes. Standardization uses pooled SDs within ROI, rating, and context pair. Confidence intervals are pointwise, not multiplicity-adjusted.\n\n');
for predictor=["pleasantness","intensity"]
    fprintf(fid,'## %s\n\n| ROI | Context pair | Raw slope [95%% CI] | Standardized beta | p | q | Status |\n',predictor);
    fprintf(fid,'| --- | --- | --- | ---: | ---: | ---: | --- |\n');
    for i=find(slopes.predictor==predictor)'
        a=slopes(i,:);
        fprintf(fid,'| %s | %s | %.4f [%.4f, %.4f] | %.4f | %.4g | %.4g | %s |\n', ...
            a.roi,a.context_pair,a.slope_raw,a.ci_lower_raw,a.ci_upper_raw,a.beta_standardized,a.p_value,a.q_value,a.model_status);
    end
    fprintf(fid,'\n');
end
output=struct('formula',formula,'slopes',slopes,'fixed_effects',fixed,'diagnostics',diagnostics, ...
    'residuals',residual_table,'predictions',predictions,'models',{fitted_models}, ...
    'convergence',{convergence_details},'covariance',{covariance_details});
save(fullfile(out,'mixed_models.mat'),'output','-v7.3');
end
function saveplot(f,stem)
exportgraphics(f,stem+".png",'Resolution',180);
exportgraphics(f,stem+".pdf",'ContentType','vector'); close(f);
end
