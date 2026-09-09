function results = OX_formal_rsa(analysis,root,varargin)
% Shared RSA estimation and global condition-correspondence permutation engine.
p=inputParser;
p.addParameter('NumPermutations',5000,@(x) isscalar(x)&&x>=1&&x==floor(x));
p.addParameter('PermutationSeed',1001,@(x) isscalar(x)&&x>=0&&x==floor(x));
p.addParameter('OutputDir','',@(x) ischar(x)||isstring(x));
p.addParameter('MakePlots',true,@islogical);
p.addParameter('RunMixedModels',true,@islogical);
p.addParameter('SubjectIDs',2:6,@(x) isnumeric(x)&&isvector(x)&&all(ismember(x,2:6))&&numel(unique(x))==numel(x));
p.addParameter('ROINames',["AON","PirF","PirT","olfAMG","olfOFC"],@(x) isstring(x)||iscellstr(x));
p.addParameter('NeuralDir',fullfile(root,'RDMs','neural'),@(x) ischar(x)||isstring(x));
p.addParameter('RequireExploratoryReference',true,@islogical);
p.parse(varargin{:}); opts=p.Results;
assert(ismember(string(analysis),["omnibus","displacement"]));
names=["omnibus_RSA","rating_displacement_RSA"];
out=string(opts.OutputDir);
if strlength(out)==0, out=fullfile(root,'RDMs',names(1+(string(analysis)=="displacement"))); end
for folder=["","subject_results","tables","figures"]
    if ~isfolder(fullfile(out,folder)), mkdir(fullfile(out,folder)); end
end
subjects=opts.SubjectIDs(:)'; rois=string(opts.ROINames(:))';
nS=numel(subjects); nR=numel(rois);
assert(nS>=2 && nR>=1 && numel(unique(rois))==nR);
contexts=["PERSON","FOOD","LOCATION","CONTROL"];
rdmroot=fullfile(root,'RDMs'); exploratory=fullfile(rdmroot,'exploratory_RSA_results','tables');
C=load(fullfile(rdmroot,'categorical_RDMs.mat'));
ui=find(triu(true(80),1)); [ii,jj]=ind2sub([80 80],ui);
pair_indices=nchoosek(1:4,2); di=zeros(120,1); dj=di; pair=strings(120,1); odor=repmat((1:20)',6,1);
for k=1:6
    rows=(k-1)*20+(1:20); di(rows)=(pair_indices(k,1)-1)*20+(1:20);
    dj(rows)=(pair_indices(k,2)-1)*20+(1:20);
    pair(rows)=contexts(pair_indices(k,1))+"_vs_"+contexts(pair_indices(k,2));
end
dispidx=sub2ind([80 80],di,dj);
pair_names=unique(pair,'stable'); pairdummy=double(pair==pair_names(2:end)');
odordummy=double(odor==(2:20));
if string(analysis)=="omnibus"
    model_names="omnibus"; predictors=["odor","context","pleasantness","intensity"];
    reference=readtable(fullfile(exploratory,'analysis1_omnibus_RSA.csv'),'TextType','string');
else
    model_names=["context_pair_plus_odor","context_pair"]; predictors=["pleasantness","intensity"];
    reference=readtable(fullfile(exploratory,'analysis6_same_odor_displacement.csv'),'TextType','string');
end
nP=numel(predictors); nB=opts.NumPermutations;
blank=struct('name',"",'predictors',predictors,'beta',nan(nS,nR,nP), ...
    'null_beta',nan(nS,nR,nP,nB),'subject_statistics',table(),'group_statistics',table(), ...
    'group_null_beta',[],'nuisance_coefficients',table(),'fit_statistics',table());
models=repmat(blank,1,numel(model_names));
for m=1:numel(models), models(m).name=model_names(m); end
permutation_mappings=zeros(80,nB,nS,'uint8');
provenance=cell(nS,1); observations=table();
for s=1:nS
    sid=subjects(s); B=load(fullfile(rdmroot,"subj_"+sid+"_behavioral_RDMs.mat"));
    N=load(fullfile(opts.NeuralDir,"subj_"+sid+"_neural_RDMs.mat"),'results'); N=N.results;
    for metadata={B.condition_metadata,N.condition_metadata,C.condition_metadata}
        M=metadata{1};
        assert(isequal(double(M.condition_index(:)),(1:80)'));
        assert(isequal(double(M.condition_odor(:)),repmat((1:20)',4,1)));
        assert(isequal(string(M.condition_context(:)),repelem(contexts(:),20)));
    end
    assert(all(ismember(C.D_odor(:),[0 1])) && all(ismember(C.D_context(:),[0 1])));
    provenance{s}=struct('subject_id',sid,'preprocessing',N.preprocessing, ...
        'roi_metadata',N.roi_metadata,'neural_input',fullfile(opts.NeuralDir,"subj_"+sid+"_neural_RDMs.mat"), ...
        'behavior_input',fullfile(rdmroot,"subj_"+sid+"_behavioral_RDMs.mat"));
    stream=RandStream('mt19937ar','Seed',opts.PermutationSeed+sid*100000);
    perms=zeros(80,nB);
    for k=1:nB, perms(:,k)=randperm(stream,80)'; end
    assert(isequal(sort(perms,1),repmat((1:80)',1,nB)));
    permutation_mappings(:,:,s)=uint8(perms);
    if s>1, assert(~isequal(permutation_mappings(:,:,s),permutation_mappings(:,:,s-1))); end
    if string(analysis)=="omnibus"
        X=zscore([C.D_odor(ui),C.D_context(ui),B.D_pleasantness(ui),B.D_intensity(ui)]);
        base=[ones(3160,1),X]; designs={base}; terms={"intercept"};
        row_i=ii; row_j=jj; observed_indices=ui;
    else
        X=zscore([B.D_pleasantness(dispidx),B.D_intensity(dispidx)]);
        base=[ones(120,1),X,pairdummy]; designs={[base,odordummy],base};
        pairterms="context_pair_"+pair_names(2:end)';
        terms={["intercept",pairterms,"odor_"+string(2:20)],["intercept",pairterms]};
        row_i=di; row_j=dj; observed_indices=dispidx;
    end
    inverse=cell(size(designs));
    for m=1:numel(models)
        assert(all(isfinite(designs{m}(:))) && rank(designs{m})==size(designs{m},2));
        [Q,R]=qr(designs{m},0); inverse{m}=R\Q';
    end
    for r=1:nR
        roi=N.roi_results(strcmp(string({N.roi_results.roi_name}),rois(r)));
        assert(isscalar(roi)); D=double(roi.D_neural_simple);
        assert(isequal(size(D),[80 80])&&all(isfinite(D(:))));
        assert(max(abs(D-D'),[],'all')<1e-10 && max(abs(diag(D)))<1e-10);
        raw=D(observed_indices); assert(std(raw)>sqrt(eps)); y=zscore(raw);
        if string(analysis)=="displacement"
            observations=[observations;table(repmat(sid,120,1),repmat(rois(r),120,1),odor,pair,raw, ...
                B.D_pleasantness(dispidx),B.D_intensity(dispidx), ...
                'VariableNames',{'subject_id','roi','odor_id','context_pair','neural_distance', ...
                'abs_delta_pleasantness','abs_delta_intensity'})]; %#ok<AGROW>
        end
        for m=1:numel(models)
            A=designs{m}; b=A\y; models(m).beta(s,r,:)=b(2:nP+1);
            refmask=reference.subject_id==sid & reference.roi==rois(r) & ...
                reference.neural_metric=="simple" & reference.version=="linear";
            if string(analysis)=="displacement", refmask=refmask & reference.model==model_names(m); end
            assert(nnz(refmask)<=1);
            reference_error=NaN;
            if any(refmask)
                ref=reference{refmask,cellstr("beta_"+predictors)};
                reference_error=max(abs(ref-b(2:nP+1)'));
                assert(reference_error<1e-10,'Exploratory beta mismatch.');
            else
                assert(~opts.RequireExploratoryReference,'Missing exploratory reference.');
            end
            independent_error=max(abs(pinv(A)*y-b));
            assert(independent_error<1e-10,'Independent SVD fit mismatch.');
            nuisance=b([1,nP+2:numel(b)]);
            models(m).nuisance_coefficients=[models(m).nuisance_coefficients; ...
                table(repmat(sid,numel(nuisance),1),repmat(rois(r),numel(nuisance),1), ...
                terms{m}(:),nuisance,'VariableNames',{'subject_id','roi','term','coefficient'})];
            sse=sum((y-A*b).^2); r2=1-sse/sum(y.^2); adj=1-(1-r2)*(numel(y)-1)/(numel(y)-numel(b));
            models(m).fit_statistics=[models(m).fit_statistics; ...
                table(sid,rois(r),numel(y),0,rank(A),r2,adj,reference_error,independent_error, ...
                'VariableNames',{'subject_id','roi','n_valid','n_dropped','design_rank','r2','adjusted_r2','exploratory_max_beta_error','independent_svd_max_beta_error'})];
        end
        for first=1:250:nB
            kk=first:min(first+249,nB);
            indices=sub2ind([80 80],perms(row_i,kk),perms(row_j,kk));
            Y=D(indices); assert(all(std(Y,0,1)>sqrt(eps))); Y=zscore(Y,0,1);
            for m=1:numel(models)
                bb=inverse{m}*Y;
                models(m).null_beta(s,r,:,kk)=reshape(bb(2:nP+1,:),[1,1,nP,numel(kk)]);
            end
        end
    end
    fprintf('[%s] subject %d: %d permutations across %d ROIs complete.\n',analysis,sid,nB,nR);
end
for m=1:numel(models)
    model=models(m); assert(all(isfinite(model.null_beta(:))) && all(isfinite(model.beta(:))));
    group_null=reshape(mean(model.null_beta,1),nR,nP,nB);
    group_beta=reshape(mean(model.beta,1),nR,nP);
    st=table(); gt=table();
    for s=1:nS
        part=table();
        for r=1:nR
            for j=1:nP
                nb=reshape(model.null_beta(s,r,j,:),[],1); obs=model.beta(s,r,j);
                part=[part;table(subjects(s),rois(r),predictors(j),obs,OX_rsa_empirical_p(obs,nb), ...
                    mean(nb),std(nb),nB,'VariableNames',{'subject_id','roi','predictor','beta','p_value','null_mean','null_sd','n_permutations'})]; %#ok<AGROW>
            end
        end
        part.q_value=OX_rsa_bh(part.p_value); st=[st;part]; %#ok<AGROW>
    end
    null_table=table();
    for r=1:nR
        for j=1:nP
            nb=reshape(group_null(r,j,:),[],1); obs=group_beta(r,j);
            assert(max(abs(nb-reshape(sum(model.null_beta(:,r,j,:),1)/nS,[],1)))<1e-12);
            gt=[gt;table(rois(r),predictors(j),obs,OX_rsa_empirical_p(obs,nb),mean(nb),std(nb),nS,nB, ...
                'VariableNames',{'roi','predictor','mean_beta','p_value','null_mean','null_sd','n_subjects','n_permutations'})]; %#ok<AGROW>
            null_table=[null_table;table((1:nB)',repmat(rois(r),nB,1),repmat(predictors(j),nB,1),nb, ...
                'VariableNames',{'permutation','roi','predictor','mean_null_beta'})]; %#ok<AGROW>
        end
    end
    gt.q_value=OX_rsa_bh(gt.p_value);
    models(m).subject_statistics=st; models(m).group_statistics=gt; models(m).group_null_beta=group_null;
    prefix=fullfile(out,'tables',model.name);
    writetable(st,prefix+"_subject_statistics.csv"); writetable(gt,prefix+"_group_statistics.csv");
    writetable(null_table,prefix+"_group_null.csv");
    writetable(model.nuisance_coefficients,prefix+"_nuisance_coefficients.csv");
    writetable(model.fit_statistics,prefix+"_fit_qc.csv");
end
results=struct('analysis',string(analysis),'options',opts,'subject_ids',subjects, ...
    'roi_names',rois,'context_names',contexts,'models',models,'provenance',{provenance}, ...
    'permutation_mappings',permutation_mappings,'permutation_seeds',opts.PermutationSeed+subjects*100000, ...
    'output_dir',out,'condition_order',C.condition_metadata,'observations',observations);
results.inference="Two-sided global neural condition-correspondence null; equal-weight mean of the selected independent subject null draws. Not a coefficient-specific conditional null or a population random-effects test.";
results.null_dimensions="models.null_beta: subject x ROI x predictor x permutation; group_null_beta: ROI x predictor x permutation";
results.scaling="Sample-SD z-scoring; neural outcome restandardized on every permutation. Categorical omnibus predictors also standardized. Displacement dummies not standardized.";
results.input_qc=readtable(fullfile(exploratory,'QC_input_status.csv'),'TextType','string');
results.input_qc=results.input_qc(ismember(results.input_qc.subject_id,subjects),:);
writetable(results.input_qc,fullfile(out,'tables','input_qc.csv'));
writetable(C.condition_metadata,fullfile(out,'tables','condition_order.csv'));
if ~isempty(observations), writetable(observations,fullfile(out,'tables','observations.csv')); end
for s=1:nS
    subject_result=struct('subject_id',subjects(s),'roi_names',rois,'permutations',permutation_mappings(:,:,s), ...
        'seed',results.permutation_seeds(s),'provenance',provenance{s});
    for m=1:numel(models)
        subject_result.models(m).name=models(m).name;
        subject_result.models(m).predictors=predictors;
        subject_result.models(m).beta=reshape(models(m).beta(s,:,:),nR,nP);
        subject_result.models(m).null_beta=reshape(models(m).null_beta(s,:,:,:),nR,nP,nB);
        subject_result.models(m).statistics=models(m).subject_statistics(models(m).subject_statistics.subject_id==subjects(s),:);
    end
    save(fullfile(out,'subject_results',"subj_"+subjects(s)+"_results.mat"),'subject_result','-v7.3');
end
save(fullfile(out,'results.mat'),'results','-v7.3');
if opts.MakePlots, OX_plot_formal_rsa(results); end
if string(analysis)=="displacement" && opts.RunMixedModels
    results.mixed_models=OX_rsa_context_displacement_lme(observations,out,opts.MakePlots);
    save(fullfile(out,'results.mat'),'results','-v7.3');
end
OX_rsa_write_readme(results);
fprintf('[%s] Complete: %s\n',analysis,out);
end
