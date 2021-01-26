function[P2,N,err_p,err_n] = isometric_nrsfm(Pgth,Ngth,qgth,I1u,I1v,I2u,I2v,J21a,J21b,J21c,J21d,J12a,J12b,J12c,J12d,H21uva,H21uvb,sm)

idx = 1:length(qgth); % the images the algorithm needs to be evaluated on (in order)
J21a(abs(J21a)<.001)= 0; J21b(abs(J21b)<.001)= 0; J21c(abs(J21c)<.001)= 0;J21d(abs(J21d)<.001)= 0;
J12a(abs(J12a)<.001)= 0; J12b(abs(J12b)<.001)= 0; J12c(abs(J12c)<.001)= 0;J12d(abs(J12d)<.001)= 0;
% coeff of k1       % coeff of k2        % constant term
T1_12_k1 = -J21c;   T1_12_k2 = -J21d;    T1_12_c = (J12a.*H21uva + J12c.*H21uvb);%(H21vvb./J21d)/2;
T2_12_k1 = -J21a;   T2_12_k2 = -J21b;    T2_12_c = (J12b.*H21uva + J12d.*H21uvb);%(H21uub./J21b)/2;

% k1b = -T2_12 = a*k1 + b*k2 + t1;
% k2b = -T1_12 = c*k1 + d*k2 + t2;
a = -T2_12_k1; b = -T2_12_k2; c = -T1_12_k1; d = -T1_12_k2; t1 = -T2_12_c; t2 = -T1_12_c;

% write metric tensor with k1b and k2b: g11, g12, g22
e = 1+ I2u.^2 + I2v.^2;
u = I2u;
v = I2v;

% Metric tensor (see equation 26 in the paper)
% pullback metric tensor (see equation 8 in the paper): g11b,g12b,g22b
e1 = repmat(1+ I1u.^2 + I1v.^2,length(idx)-1,1);
u1 = repmat(I1u ,length(idx)-1,1);
v1 = repmat(I1v ,length(idx)-1,1);

% Reconstruction equations (see equation 29  in paper)
% eq1 = g11*g22b-g11b*g22) eq2 = g12*g22b-g12b*g22

% Minimise sum of squares expression 
% eq = eq1^2 + eq2^2

eq = create_polynomial_coefficients(a,b,c,d,t1,t2,e,e1,u,u1,v,v1);
res = solve_polynomial(eq);

% recover first order derivatives on rest of the surfaces
k1_all = [res(:,1)';a.*repmat(res(:,1)',length(idx)-1,1) + b.*repmat(res(:,2)',length(idx)-1,1) + t1];
k2_all = [res(:,2)';c.*repmat(res(:,1)',length(idx)-1,1) + d.*repmat(res(:,2)',length(idx)-1,1) + t2];

u_all = [I1u;I2u]; v_all = [I1v;I2v];

% find normals on all surfaces N= [N1;N2;N3]
N1 = k1_all; N2 = k2_all; N3 = 1-u_all.*k1_all-v_all.*k2_all;
n = sqrt(N1.^2+N2.^2+N3.^2);
N1 = N1./n ; N2 = N2./n; N3 = N3./n;

N = [N1(:),N2(:),N3(:)]';
N_res = reshape(N(:),3*length(idx),length(u_all));


% find indices with no solution
idx = find(res(:,1)==0);
N_res(:,idx) = []; u_all(:,idx) = []; v_all(:,idx) = [];


% Integrate normals to find depth
 P_grid=calculate_depth(N_res,u_all,v_all,sm);

% compare with ground truth
[P2,err_p] = compare_with_Pgth(P_grid,u_all,v_all,qgth,Pgth);
[N,err_n] = compare_with_Ngth(P2,qgth,Ngth);

%plot results
% for i=1:size(u_all,1)
%      figure(i)
%     %plot3(Pgth(3*(i-1)+1,:),Pgth(3*(i-1)+2,:),Pgth(3*(i-1)+3,:),'go');
%     hold on;
%      % plot3(P2(3*(i-1)+1,:),P2(3*(i-1)+2,:),P2(3*(i-1)+3,:),'ro');
% %       plot3(P_grid(3*(i-1)+1,:),P_grid(3*(i-1)+2,:),P_grid(3*(i-1)+3,:),'bo');
%       quiver3(P_grid(3*(i-1)+1,:),P_grid(3*(i-1)+2,:),P_grid(3*(i-1)+3,:),N_res(3*(i-1)+1,:),N_res(3*(i-1)+2,:),N_res(3*(i-1)+3,:));
%       %quiver3(Pgth(3*(i-1)+1,:),Pgth(3*(i-1)+2,:),Pgth(3*(i-1)+3,:),Ngth(3*(i-1)+1,:),Ngth(3*(i-1)+2,:),Ngth(3*(i-1)+3,:));
%     hold off;
%     axis equal;
% end
% mean(err_p')
% mean(err_n')
