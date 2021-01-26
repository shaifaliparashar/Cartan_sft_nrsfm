function [J21,Huu21,Huv21,Hvv21,qw2] = create_warps(P2,q,n,par)

er = 1e-4;
nC = 100;
t= 1e-3;

% 1:  Eta_21 derivatives
for i=2:n
    i
    umin = min(q{i}(1,:))-t; umax = max(q{i}(1,:))+t;
    vmin = min(q{i}(2,:))-t; vmax = max(q{i}(2,:))+t;
    bbs = bbs_create(umin, umax, nC, vmin, vmax, nC, 2);
    
    coloc = bbs_coloc(bbs, q{i}(1,:), q{i}(2,:));
    lambdas = er*ones(nC-3, nC-3);
    bending = bbs_bending(bbs, lambdas);
    
    % get control points for i to j warp
    cpts = (coloc'*coloc + bending) \ (coloc'*q{1}(1:2,:)');
    ctrlpts = cpts';
            [xv,yv]=meshgrid(linspace(umin,umax,100),linspace(vmin,vmax,100));
            ctrlpts = optimPanalSchwarz(bbs,ctrlpts,q{i}',q{1}',[xv(:),yv(:)],par(i));
            ctrlpts=ctrlpts';
    qw2{i-1} = bbs_eval(bbs,ctrlpts,q{i}(1,:)',q{i}(2,:)',0,0);
    error=sqrt(mean((qw2{i-1}(1,:)-q{1}(1,:)).^2+(qw2{i-1}(2,:)-q{1}(2,:)).^2));
    dqu = bbs_eval(bbs, ctrlpts, q{i}(1,:)',q{i}(2,:)',1,0);
    dqv = bbs_eval(bbs, ctrlpts, q{i}(1,:)',q{i}(2,:)',0,1);
    dquv = bbs_eval(bbs,ctrlpts,q{i}(1,:)',q{i}(2,:)',1,1);
    dquu = bbs_eval(bbs, ctrlpts, q{i}(1,:)',q{i}(2,:)',2,0);
    dqvv = bbs_eval(bbs, ctrlpts, q{i}(1,:)',q{i}(2,:)',0,2);
    J21{i-1} = [dqu; dqv];
    Huu21{i-1} = dquu;
    Huv21{i-1} = dquv;
    Hvv21{i-1} = dqvv;
    disp([sprintf('[ETA] Internal Rep error = %f',error)]);
    %Visualize Point Registration Error
%     [xv,yv]=meshgrid(linspace(umin+t,umax-t,20),linspace(vmin+t,vmax-t,20));
%     qv = bbs_eval(bbs,ctrlpts,xv(:),yv(:),0,0);
%     figure;
%     plot(q{1}(1,:),q{1}(2,:),'ro');
%     hold on;
%     plot(qw2{i-1}(1,:),qw2{i-1}(2,:),'b*');
%     mesh(reshape(qv(1,:),size(xv)),reshape(qv(2,:),size(xv)),zeros(size(xv)));
%     axis equal
%     hold off;
end

% % 2:  Eta_12 derivatives
% for i=2:n
%     
%     umin = min(q{1}(1,:))-t; umax = max(q{1}(1,:))+t;
%     vmin = min(q{1}(2,:))-t; vmax = max(q{1}(2,:))+t;
%     bbs = bbs_create(umin, umax, nC, vmin, vmax, nC, 2);
%     
%     coloc = bbs_coloc(bbs, q{1}(1,:), q{1}(2,:));
%     lambdas = er*ones(nC-3, nC-3);
%     bending = bbs_bending(bbs, lambdas);
%     
%     % get control points for i to j warp
%     cpts = (coloc'*coloc + bending) \ (coloc'*q{i}(1:2,:)');
%     ctrlpts = cpts';
% %             [xv,yv]=meshgrid(linspace(umin,umax,100),linspace(vmin,vmax,100));
% %             ctrlpts = optimPanalSchwarz(bbs,ctrlpts,q{1}',q{i}',[xv(:),yv(:)],par(i));
% %             ctrlpts=ctrlpts';
%     qw{i-1} = bbs_eval(bbs,ctrlpts,q{1}(1,:)',q{1}(2,:)',0,0);
%     error=sqrt(mean((qw{i-1}(1,:)-q{i}(1,:)).^2+(qw{i-1}(2,:)-q{i}(2,:)).^2));
%     dqu = bbs_eval(bbs, ctrlpts, q{1}(1,:)',q{1}(2,:)',1,0);
%     dqv = bbs_eval(bbs, ctrlpts, q{1}(1,:)',q{1}(2,:)',0,1);
%     J12{i-1} = [dqu; dqv];
%     disp([sprintf('[ETA] Internal Rep error = %f',error)]);
%     %Visualize Point Registration Error
%     [xv,yv]=meshgrid(linspace(umin+t,umax-t,20),linspace(vmin+t,vmax-t,20));
%     qv = bbs_eval(bbs,ctrlpts,xv(:),yv(:),0,0);
%     figure;
%     plot(q{i}(1,:),q{i}(2,:),'ro');
%     hold on;
%     plot(qw{i-1}(1,:),qw{i-1}(2,:),'b*');
%     mesh(reshape(qv(1,:),size(xv)),reshape(qv(2,:),size(xv)),zeros(size(xv)));
%     axis equal
%     hold off;
% end
% 
