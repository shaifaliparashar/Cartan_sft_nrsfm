          function [Et,Ft,Gt,x1,x2,N,Pest]= create_ground_truth_normals_and_metric_tensor(P,q,u,v,ub,vb)

% create Ground truth
er = 1e-5;
t= 1e-3;
nC = 100;

umin = min(min(q{1}(1,:)),min(u))-t; umax = max(max(q{1}(1,:)),max(u))+t;
vmin = min(min(q{1}(2,:)),min(v))-t; vmax = max(max(q{1}(2,:)),max(v))+t;
bbs = bbs_create(umin, umax, nC, vmin, vmax, nC, 3);
coloc = bbs_coloc(bbs, q{1}(1,:), q{1}(2,:));
lambdas = er*ones(nC-3, nC-3);
bending = bbs_bending(bbs, lambdas);
cpts = (coloc'*coloc + bending) \ (coloc'*P(1:3,:)');
ctrlpts = cpts';
Pest{1} = bbs_eval(bbs, ctrlpts, u',v',0,0);
dqu = bbs_eval(bbs, ctrlpts, u',v',1,0);
dqu = dqu./repmat(sqrt(sum(dqu.^2)),3,1);
dqv = bbs_eval(bbs, ctrlpts, u',v',0,1);
dqv = dqv./repmat(sqrt(sum(dqv.^2)),3,1);
x1 = -dqu(3,:);
x2 = -dqv(3,:);
N{1} = cross(dqu,dqv);
N{1} = N{1}./repmat(sqrt(sum(N{1}.^2)),3,1);
Et = (1+ u.^2 + v.^2).*(x1.^2) + 1 - 2.*u.*x1;
Ft = (1+ u.^2 + v.^2).*(x2.^2) + 1 - 2.*v.*x2;
Gt = (1+ u.^2 + v.^2).*(x1.*x2) - u.*x2 - v.*x1;

for i = 2:length(q)
    umin = min(min(q{i}(1,:)),min(ub(i-1,:)))-t; umax = max(max(q{i}(1,:)),max(ub(i-1,:)))+t;
    vmin = min(min(q{i}(2,:)),min(vb(i-1,:)))-t; vmax = max(max(q{i}(2,:)),max(vb(i-1,:)))+t;
    bbs = bbs_create(umin, umax, nC, vmin, vmax, nC, 3);
    coloc = bbs_coloc(bbs, q{i}(1,:), q{i}(2,:));
    lambdas = er*ones(nC-3, nC-3);
    bending = bbs_bending(bbs, lambdas);
    cpts = (coloc'*coloc + bending) \ (coloc'*P(3*(i-1) + 1:3*(i-1)+3,:)');
    ctrlpts = cpts';
    Pest{i} = bbs_eval(bbs, ctrlpts, ub(i-1,:)',vb(i-1,:)',0,0);
    dqu = bbs_eval(bbs, ctrlpts, ub(i-1,:)',vb(i-1,:)',1,0);
    dqu = dqu./repmat(sqrt(sum(dqu.^2)),3,1);
    dqv = bbs_eval(bbs, ctrlpts, ub(i-1,:)',vb(i-1,:)',0,1);
    dqv = dqv./repmat(sqrt(sum(dqv.^2)),3,1);
    N{i} = cross(dqu,dqv);
    N{i} = N{i}./repmat(sqrt(sum(N{i}.^2)),3,1);
end








