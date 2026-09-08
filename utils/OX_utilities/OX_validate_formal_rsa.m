function report = OX_validate_formal_rsa(root)
% Validate saved production artifacts against direct refits and CSV exports.
if nargin==0, root=fileparts(fileparts(fileparts(mfilename('fullpath')))); end
O=load(fullfile(root,'RDMs','omnibus_RSA','results.mat')); A=O.results;
O=load(fullfile(root,'RDMs','rating_displacement_RSA','results.mat')); B=O.results;
assert(isequal(A.permutation_mappings,B.permutation_mappings));
assert(isequal(A.subject_ids,2:6)&&isequal(B.subject_ids,2:6));
assert(A.options.NumPermutations==5000 && B.options.NumPermutations==5000);
assert(all(all(sort(double(A.permutation_mappings),1)==repmat((1:80)',1,5000,5)),'all'));
C=load(fullfile(root,'RDMs','categorical_RDMs.mat')); ui=find(triu(true(80),1));
pairdummy=kron(eye(6),ones(20,1)); pairdummy=pairdummy(:,2:end);
odordummy=repmat(eye(20),6,1); odordummy=odordummy(:,2:end);
cp=nchoosek(1:4,2); di=[]; dj=[];
for k=1:6
    di=[di;(cp(k,1)-1)*20+(1:20)']; dj=[dj;(cp(k,2)-1)*20+(1:20)']; %#ok<AGROW>
end
idx=sub2ind([80 80],di,dj); max_error=0; n_refits=0;
for s=1:5
    sid=A.subject_ids(s); N=load(fullfile(root,'RDMs','neural',"subj_"+sid+"_neural_RDMs.mat"));
    P=load(fullfile(root,'RDMs',"subj_"+sid+"_behavioral_RDMs.mat"));
    X1=[ones(3160,1),zscore([C.D_odor(ui),C.D_context(ui),P.D_pleasantness(ui),P.D_intensity(ui)])];
    X2=[ones(120,1),zscore([P.D_pleasantness(idx),P.D_intensity(idx)]),pairdummy];
    X3=[X2,odordummy];
    for r=1:5
        R=N.results.roi_results(strcmp(string({N.results.roi_results.roi_name}),A.roi_names(r)));
        for k=[1 2 19 5000]
            perm=double(A.permutation_mappings(:,k,s)); D=R.D_neural_simple(perm,perm);
            b1=X1\zscore(D(ui)); b2=X2\zscore(D(idx)); b3=X3\zscore(D(idx));
            errors=[b1(2:5)-reshape(A.models(1).null_beta(s,r,:,k),[],1); ...
                b2(2:3)-reshape(B.models(2).null_beta(s,r,:,k),[],1); ...
                b3(2:3)-reshape(B.models(1).null_beta(s,r,:,k),[],1)];
            max_error=max(max_error,max(abs(errors))); n_refits=n_refits+3;
        end
    end
end
assert(max_error<1e-10);
for item={A,B}
    R=item{1};
    for m=1:numel(R.models)
        M=R.models(m); expected=reshape(mean(M.null_beta,1),size(M.group_null_beta));
        assert(max(abs(expected-M.group_null_beta),[],'all')<1e-12);
        for s=1:5
            S=M.subject_statistics(M.subject_statistics.subject_id==R.subject_ids(s),:);
            for i=1:height(S)
                r=find(R.roi_names==S.roi(i)); j=find(M.predictors==S.predictor(i));
                assert(abs(S.p_value(i)-OX_rsa_empirical_p(S.beta(i),M.null_beta(s,r,j,:)))<1e-12);
            end
            assert(max(abs(S.q_value-OX_rsa_bh(S.p_value)))<1e-12);
        end
        G=readtable(fullfile(R.output_dir,'tables',M.name+"_group_statistics.csv"),'TextType','string');
        for i=1:height(G)
            r=find(R.roi_names==G.roi(i)); j=find(M.predictors==G.predictor(i));
            assert(abs(G.mean_beta(i)-mean(M.beta(:,r,j)))<1e-12);
            assert(abs(G.p_value(i)-OX_rsa_empirical_p(G.mean_beta(i),M.group_null_beta(r,j,:)))<1e-12);
        end
        assert(max(abs(G.q_value-OX_rsa_bh(G.p_value)))<1e-12);
    end
end
assert(isfield(B,'mixed_models') && height(B.mixed_models.slopes)==30);
assert(B.mixed_models.grouping=="semantic_context_vs_rest");
assert(all(B.mixed_models.slopes.n_observations==300) && all(B.mixed_models.slopes.n_context_pairs==3));
assert(height(B.mixed_models.membership)==9);
for c=B.mixed_models.contexts
    pairs=B.mixed_models.membership.context_pair(B.mixed_models.membership.focal_context==c);
    assert(numel(pairs)==3 && all(any(split(pairs,'_vs_')==c,2)));
    assert(nnz(contains(pairs,'CONTROL'))==1);
end
assert(all(isfinite(B.mixed_models.slopes.p_value)));
assert(max(abs(B.mixed_models.slopes.q_value-OX_rsa_bh(B.mixed_models.slopes.p_value)))<1e-12);
report=table(["Shared 5000 bijections, independent subject seeds";"Direct refits of 300 permuted models"; ...
    "Individual and group empirical p-values and BH-FDR";"Equal-weight group null construction";"30 context-versus-rest mixed-model slopes and subset membership"], ...
    repmat("pass",5,1),[0;max_error;0;0;0],'VariableNames',{'check','status','max_numeric_error'});
for out=[A.output_dir,B.output_dir], writetable(report,fullfile(out,'tables','validation.csv')); end
fprintf('Production validation passed: %d direct permutation refits, max error %.3g.\n',n_refits,max_error);
end
