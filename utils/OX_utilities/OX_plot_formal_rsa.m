function OX_plot_formal_rsa(results)
% Boxcharts summarize subjects, never permutations or RDM entries.
rois=results.roi_names; subjects=results.subject_ids; nS=numel(subjects); nR=numel(rois); subject_colors=lines(5); sc=subject_colors(subjects-1,:);
figure_note = "Group signed-beta tests: * q < 0.05, ** q < 0.01, *** q < 0.001";
if isfield(results,'figure_note'), figure_note=string(results.figure_note); end
for m=1:numel(results.models)
    M=results.models(m); pred=M.predictors; pc=lines(4);
    if results.analysis=="displacement", pc=pc(3:4,:); end
    if results.analysis=="omnibus"
        f=figure('Visible','off','Color','w','Position',[50 50 1300 850]);
        t=tiledlayout(f,2,2,'TileSpacing','compact');
        for j=1:numel(pred)
            ax=nexttile(t); hold(ax,'on');
            panelmax=max(abs(M.beta(:,:,j)),[],'all','omitnan');
            for r=1:nR
                b=abs(M.beta(:,r,j)); valid=isfinite(b); offsets=linspace(-.12,.12,nS)';
                boxchart(ax,repmat(r,nnz(valid),1),b(valid),'BoxFaceColor',pc(j,:),'MarkerStyle','none');
                scatter(ax,r+offsets(valid),b(valid),40,sc(valid,:),'filled');
                add_q_star(ax,r,b,M,rois(r),pred(j),panelmax);
            end
            ylim(ax,[0,max(panelmax,eps)*1.22]);
            ax.XTick=1:nR; ax.XTickLabel=rois+" (n="+string(sum(isfinite(M.beta(:,:,j)),1))+")"; ax.TickLabelInterpreter='none'; yline(ax,0,':k');
            title(ax,pred(j)); ylabel(ax,'Absolute standardized beta');
        end
        h=gobjects(numel(subjects),1);
        for s=1:numel(subjects), h(s)=scatter(ax,nan,nan,40,sc(s,:),'filled'); end
        lg=legend(ax,h,"Subject "+string(subjects),'Orientation','horizontal'); lg.Layout.Tile='south';
        title(t,["Omnibus RSA | simple distance"; figure_note],'Interpreter','none');
        saveplot(f,fullfile(results.output_dir,'figures','betas_by_predictor'));
        f=figure('Visible','off','Color','w','Position',[50 50 1450 850]);
        ncols=min(3,nR);
        t=tiledlayout(f,ceil(nR/ncols),ncols,'TileSpacing','compact');
        for r=1:nR
            ax=nexttile(t); hold(ax,'on');
            panelmax=max(abs(M.beta(:,r,:)),[],'all','omitnan');
            for j=1:numel(pred)
                b=abs(M.beta(:,r,j)); valid=isfinite(b); offsets=linspace(-.12,.12,nS)';
                boxchart(ax,repmat(j,nnz(valid),1),b(valid),'BoxFaceColor',pc(j,:),'MarkerStyle','none');
                scatter(ax,j+offsets(valid),b(valid),40,sc(valid,:),'filled');
                add_q_star(ax,j,b,M,rois(r),pred(j),panelmax);
            end
            ylim(ax,[0,max(panelmax,eps)*1.22]);
            ax.XTick=1:numel(pred); ax.XTickLabel=pred; xtickangle(ax,20);
            yline(ax,0,':k'); title(ax,rois(r)+" (n="+string(nnz(isfinite(M.beta(:,r,1))))+")",'Interpreter','none'); ylabel(ax,'Absolute standardized beta');
        end
        h=gobjects(nS,1);
        for s=1:nS, h(s)=scatter(ax,nan,nan,40,sc(s,:),'filled'); end
        lg=legend(ax,h,"Subject "+string(subjects),'Orientation','horizontal'); lg.Layout.Tile='south';
        title(t,["Omnibus RSA | simple distance"; figure_note],'Interpreter','none');
        saveplot(f,fullfile(results.output_dir,'figures','betas_by_roi'));
    else
        f=figure('Visible','off','Color','w','Position',[50 50 1350 600]);
        ax=axes(f); hold(ax,'on');
        for r=1:nR
            for j=1:2
                pos=r+(j-1.5)*.34; b=M.beta(:,r,j);
                boxchart(ax,repmat(pos,nS,1),b,'BoxFaceColor',pc(j,:),'BoxWidth',.27,'MarkerStyle','none','HandleVisibility','off');
                scatter(ax,pos+linspace(-.075,.075,nS)',b,38,sc,'filled','HandleVisibility','off');
            end
        end
        handles=gobjects(nS+2,1); labels=[pred,"Subject "+string(subjects)];
        for j=1:2, handles(j)=patch(ax,nan,nan,pc(j,:)); end
        for s=1:numel(subjects), handles(s+2)=scatter(ax,nan,nan,35,sc(s,:),'filled'); end
        legend(ax,handles,labels,'Location','eastoutside','AutoUpdate','off');
        ax.XTick=1:nR; ax.XTickLabel=rois; ax.TickLabelInterpreter='none'; xlim(ax,[.5 nR+.5]);
        yline(ax,0,':k'); ylabel(ax,'Standardized beta');
        title(ax,["Rating displacement | simple distance";replace(M.name,'_',' ') + " adjustment"]);
        saveplot(f,fullfile(results.output_dir,'figures',"coefficients_"+M.name));
    end
end
end
function add_q_star(ax,x,b,M,roi,predictor,panelmax)
G=M.group_statistics; idx=G.roi==roi & G.predictor==predictor;
assert(nnz(idx)==1 && isfinite(G.q_value(idx)),'Missing group q-value.');
q=G.q_value(idx); label="";
if q<.001, label="***"; elseif q<.01, label="**"; elseif q<.05, label="*"; end
if strlength(label)>0
    text(ax,x,max(b,[],'omitnan')+.04*max(panelmax,eps),label,'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom','FontSize',14,'FontWeight','bold','Interpreter','none');
end
end
function subjectkey(ax,subjects,colors)
hold(ax,'on'); h=gobjects(numel(subjects),1);
for s=1:numel(subjects), h(s)=scatter(ax,nan,nan,40,colors(s,:),'filled'); end
legend(ax,h,"Subject "+string(subjects),'Location','northwest'); axis(ax,'off');
end
function saveplot(f,stem)
exportgraphics(f,stem+".png",'Resolution',180);
exportgraphics(f,stem+".pdf",'ContentType','vector'); close(f);
end
