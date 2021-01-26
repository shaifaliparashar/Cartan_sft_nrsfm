function [N,err_n] = compare_with_Ngth(P,q,Ng)
er = 1e-4;
nC = 40;
for i = 1: size(P,1)/3
    if iscell(q)
        q1 = q{1,i};
    else
        q1 =q;
    end
    umin=min(q1(1,:))-0.1;umax=max(q1(1,:))+0.1;
    vmin=min(q1(2,:))-0.1;vmax=max(q1(2,:))+0.1;
    bbs = bbs_create(umin, umax, nC, vmin, vmax, nC, 3);
    coloc = bbs_coloc(bbs, q1(1,:), q1(2,:));
    lambdas = er*ones(nC-3, nC-3);
    bending = bbs_bending(bbs, lambdas);
    cpts = (coloc'*coloc + bending) \ (coloc'*P(3*(i-1)+1:3*(i-1)+3,:)');
    ctrlpts = cpts';
    qw = bbs_eval(bbs, ctrlpts, q1(1,:)',q1(2,:)',0,0);
    error = sqrt(mean(sum((qw-P(3*(i-1)+1:3*(i-1)+3,:)).^2)));
    dqu = bbs_eval(bbs, ctrlpts, q1(1,:)',q1(2,:)',1,0);
    dqv = bbs_eval(bbs, ctrlpts, q1(1,:)',q1(2,:)',0,1);
    N(3*(i-1)+1:3*(i-1)+3,:) = -cross(dqu,dqv);
    N(3*(i-1)+1:3*(i-1)+3,:) = N(3*(i-1)+1:3*(i-1)+3,:)./repmat(sqrt(sum(N(3*(i-1)+1:3*(i-1)+3,:).^2)),3,1);
    N_res = N(3*(i-1)+1:3*(i-1)+3,:);
    err_n(i,:) = acosd(dot(N_res,Ng(3*(i-1)+1:3*(i-1)+3,:)));
    
end