function OX_test_formal_rsa()
% Focused inference and matrix-index checks; no data files modified.
assert(abs(OX_rsa_empirical_p(10,1:9)-.2)<eps);
assert(abs(OX_rsa_empirical_p(0,1:9)-.2)<eps);
assert(OX_rsa_empirical_p(5,1:9)==1);
assert(OX_rsa_empirical_p(1,ones(9,1))==1);
assert(max(abs(OX_rsa_bh([.01;.04;.03;.2])-[.04;.0533333333333333;.0533333333333333;.2]))<1e-12);
assert(isequal(OX_rsa_bh([0;1]),[0;1]));
stream=RandStream('mt19937ar','Seed',201001); a=randperm(stream,80);
stream=RandStream('mt19937ar','Seed',201001); b=randperm(stream,80); assert(isequal(a,b));
stream=RandStream('mt19937ar','Seed',301001); assert(~isequal(a,randperm(stream,80)));
D=reshape(1:6400,80,80); D=(D+D')/2; D(1:81:end)=0;
P=D(a,a); assert(isequal(P,P')&&all(diag(P)==0));
[i,j]=find(triu(true(80),1)); ind=sub2ind([80 80],a(i),a(j));
assert(isequal(D(ind(:)),P(triu(true(80),1))));
for c1=1:3
    for c2=c1+1:4
        for odor=1:20
            i=(c1-1)*20+odor; j=(c2-1)*20+odor;
            assert(P(i,j)==D(a(i),a(j)));
        end
    end
end
fprintf('Formal RSA inference and permutation-index tests passed.\n');
end
