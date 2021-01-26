function [P2,Ni,err_p,err_n,P2c,Nc,err_pc,err_nc] = solve_iso_con_sft_final(Pgth,Ngth,qgth,...
    I1u,I1v,I2u,I2v,J21a,J21b,J21c,J21d,...
    H21uua,H21uub,H21vva,H21vvb,th)
% calculate ground truth normals and metric tensor % x1 and x2
[Et,Ft,Gt,x1,x2,N,Pest]= create_ground_truth_normals_and_metric_tensor(Pgth,qgth,I1u,I1v,I2u,I2v);

tt = size(J21a,1);
% create E F G
E = (J21a.^2).*repmat(Et,tt,1) + 2.*J21a.*J21b.*repmat(Gt,tt,1) + (J21b.^2).*repmat(Ft,tt,1);
G = (J21c.^2).*repmat(Et,tt,1) + 2.*J21c.*J21d.*repmat(Gt,tt,1) + (J21d.^2).*repmat(Ft,tt,1);
F = (J21a.*J21c).*repmat(Et,tt,1) + (J21a.*J21d + J21b.*J21c).*repmat(Gt,tt,1) + (J21b.*J21d).*repmat(Ft,tt,1);

A1 = J21a.*repmat(Et,tt,1) + J21b.*repmat(Gt,tt,1);
A2 = J21a.*repmat(Gt,tt,1) + J21b.*repmat(Ft,tt,1);
C1 = J21c.*repmat(Et,tt,1) + J21d.*repmat(Gt,tt,1);
C2 = J21c.*repmat(Gt,tt,1) + J21d.*repmat(Ft,tt,1);

k1 = J21a.*repmat(x1,tt,1) + J21b.*repmat(x2,tt,1)...
    - ((A1-C1).*H21uua + (A2-C2).*H21uub)./(2*(E-F));
k2 = J21c.*repmat(x1,tt,1) + J21d.*repmat(x2,tt,1)...
    - ((A1-C1).*H21vva + (A2-C2).*H21vvb)./(2*(F-G));


% calculate lambda
H1 = G.*(k1.^2)-E.*(k2.^2);
H2 = F.*(k1.^2)-E.*k1.*k2;
lam = (F.*H1-(G-E).*H2)./(((2.*I2u.*k1.*F-I2u.*k2.*E-I2v.*k1.*E).*H1)-(2.*(-I2v.*k2.*E+I2u.*k1.*G).*H2));
% D = 4.*(I2v.*k2.*E-I2u.*k1.*G).^2 -4.*eps.*(G.*(k1.^2)-E.*(k2.^2)).*(G-E);
% D(D<=0)= 0;
% al1 = (-2.*(I2v.*k2.*E-I2u.*k1.*G) + sqrt(D))./(2.*(G.*(k1.^2)-E.*(k2.^2)));
% al2 = (-2.*(I2v.*k2.*E-I2u.*k1.*G) - sqrt(D))./(2.*(G.*(k1.^2)-E.*(k2.^2)));
% al1(D==0)= 0;
% al2(D==0)= 0;
% Ep = eps.*(k1.^2) + 1 - 2.*I2u.*k1;
% Fp = eps.*k1.*k2 - I2u.*k2 - I1v.*k1;
% Gp = eps.*(k2.^2) + 1 - 2.*I2v.*k2;
% b = (Ep1+Fp1+Gp1)./(E+F+G);
% Ep1 = al1.^2.*eps.*k1.^2 + 1 - 2.*I2u.*al1.*k1;
% Fp1 = al1.^2.*eps.*k1.*k2 - I2u.*al1.*k2 - I1v.*al1.*k1;
% Gp1 = al1.^2.*eps.*k2.^2 + 1 - 2.*I2v.*al1.*k2;
% b1 = (Ep1+Fp1+Gp1)./(al1.*(E+F+G));
% Ep2 = al2.^2.*eps.*k1.^2 + 1 - 2.*I2u.*al2.*k1;
% Fp2 = al2.^2.*eps.*k1.*k2 - I2u.*al2.*k2 - I1v.*al2.*k1;
% Gp2 = al2.^2.*eps.*k2.^2 + 1 - 2.*I2v.*al2.*k2;
% b2 = (Ep2+Fp2+Gp2)./(al2.*(E+F+G));

% con1a = F.*(al1.^2.*eps.*(k1.^2)+1-2.*I2u.*al1.*k1)-E.*(al1.^2.*eps.*k1.*k2-I2u.*al1.*k2-I2v.*al1.*k1);
% con1b = F.*(al1.^2.*eps.*(k2.^2)+1-2.*I2v.*al1.*k2)-G.*(al1.^2.*eps.*k1.*k2-I2u.*al1.*k2-I2v.*al1.*k1);
% con1 = abs(con1a)+abs(con1b);
% 
% con2a = F.*(al2.^2.*eps.*(k1.^2)+1-2.*I2u.*al2.*k1)-E.*(al2.^2.*eps.*k1.*k2-I2u.*al2.*k2-I2v.*al2.*k1);
% con2b = F.*(al2.^2.*eps.*(k2.^2)+1-2.*I2v.*al2.*k2)-G.*(al2.^2.*eps.*k1.*k2-I2u.*al2.*k2-I2v.*al2.*k1);
% con2 = abs(con2a)+abs(con2b);
% 
% al = al1;
% al(con1>=con2)= al2(con1>=con2);
 lam(lam==0)=1;
% solve the equations

for i = 1: size(I2u,1) % for each input image
    Nest{i} = [k1(i,:);k2(i,:);1-I2u(i,:).*k1(i,:)-I2v(i,:).*k2(i,:)];
    Nest{i} = Nest{i}./repmat(sqrt(sum(Nest{i}.^2)),3,1);
    
    Nestc{i} = [k1(i,:);k2(i,:);lam(i,:).^(-1)-I2u(i,:).*k1(i,:)-I2v(i,:).*k2(i,:)];
    Nestc{i} = Nestc{i}./repmat(sqrt(sum(Nestc{i}.^2)),3,1);
    
    u_all = I2u(i,:); v_all = I2v(i,:);

    % Integrate normals to find depth
    P_grid=calculate_depth(Nest{i},u_all,v_all,1e-2);
    P_gridc=calculate_depth(Nestc{i},u_all,v_all,1e-2);
    
    % compare with ground truth
    [P2{i},err_p{i}] = compare_with_Pgth(P_grid,u_all,v_all,qgth{i+1},Pgth(3*i+1:3*i+3,:));
    [Ni{i},err_n{i}] = compare_with_Ngth(P2{i},qgth{i+1},Ngth(3*i+1:3*i+3,:));
    
    [P2c{i},err_pc{i}] = compare_with_Pgth(P_gridc,u_all,v_all,qgth{i+1},Pgth(3*i+1:3*i+3,:));
    [Nc{i},err_nc{i}] = compare_with_Ngth(P2c{i},qgth{i+1},Ngth(3*i+1:3*i+3,:));
end
