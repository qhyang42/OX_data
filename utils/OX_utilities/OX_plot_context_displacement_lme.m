function OX_plot_context_displacement_lme(O,out)
% Six main figures: three focal contexts times two univariate ratings.
colors=lines(3);
for c=O.contexts
    pairs=O.membership.context_pair(O.membership.focal_context==c);
    for predictor=O.predictors
        f=figure('Visible','off','Color','w','Position',[50 50 1450 900]);
        t=tiledlayout(f,2,3,'TileSpacing','compact');
        d=figure('Visible','off','Color','w','Position',[50 50 1150 1500]);
        dt=tiledlayout(d,5,2,'TileSpacing','compact');
        for r=1:5
            mask=O.residuals.focal_context==c & O.residuals.predictor==predictor & O.residuals.roi==O.roi_names(r);
            S=O.residuals(mask,:); P=O.predictions(O.predictions.focal_context==c & O.predictions.predictor==predictor & O.predictions.roi==O.roi_names(r),:);
            A=O.slopes(O.slopes.focal_context==c & O.slopes.predictor==predictor & O.slopes.roi==O.roi_names(r),:);
            ax=nexttile(t,r); hold(ax,'on');
            fill(ax,[P.absolute_rating_change;flipud(P.absolute_rating_change)], ...
                [P.ci_lower;flipud(P.ci_upper)],[.3 .3 .3],'FaceAlpha',.15,'EdgeColor','none');
            for k=1:3
                idx=S.context_pair==pairs(k);
                scatter(ax,S.absolute_rating_change(idx),S.neural_distance(idx),12,colors(k,:),'filled','MarkerFaceAlpha',.35);
            end
            plot(ax,P.absolute_rating_change,P.prediction,'k-','LineWidth',2);
            xlabel(ax,"Absolute change in "+predictor); ylabel(ax,'Neural distance (1-r)');
            label=O.roi_names(r);
            if contains(A.model_status,"fit_warning"), label=label+" (covariance warning)";
            elseif contains(A.model_status,"boundary"), label=label+" (boundary variance)"; end
            title(ax,[label;sprintf('beta = %.3f, p = %.3g, q = %.3g',A.beta_standardized,A.p_value,A.q_value)]);
            da=nexttile(dt,2*r-1); scatter(da,S.conditional_fitted,S.residual,8,'filled'); yline(da,0,':k');
            title(da,O.roi_names(r)); xlabel(da,'Conditional fitted distance'); ylabel(da,'Residual');
            da=nexttile(dt,2*r); qqplot(da,S.residual); title(da,O.roi_names(r)+" residual Q-Q");
        end
        ax=nexttile(t,6); hold(ax,'on'); h=gobjects(4,1);
        for k=1:3, h(k)=scatter(ax,nan,nan,25,colors(k,:),'filled'); end
        h(4)=plot(ax,nan,nan,'k-','LineWidth',2);
        legend(ax,h,[replace(pairs,'_',' ');"Common slope; 95% confidence band"],'Location','northwest'); axis(ax,'off');
        title(t,[c+" versus other contexts | displacement vs "+predictor; ...
            "One common rating slope; pair-specific intercepts; line averages the three pair intercepts"]);
        saveplot(f,fullfile(out,'figures',"displacement_vs_"+predictor+"_"+lower(c)+"_vs_rest_mixed_model"));
        title(dt,c+" versus rest | "+predictor+" residual diagnostics");
        saveplot(d,fullfile(out,'figures',"mixed_model_diagnostics_"+predictor+"_"+lower(c)+"_vs_rest"));
    end
end
end
function saveplot(f,stem)
exportgraphics(f,stem+".png",'Resolution',180);
exportgraphics(f,stem+".pdf",'ContentType','vector'); close(f);
end
