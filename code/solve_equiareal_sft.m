 function [P2,N,err_p,err_n] = solve_equiareal_sft(Pgth,Ngth,qgth,I1u,I1v,I2u,I2v,J21a,J21b,J21c,J21d,J12a,J12b,J12c,J12d,H21uua,H21uub,H21uva,H21uvb,H21vva,H21vvb,th)
J21a(abs(J21a)<.001)= 0; J21b(abs(J21b)<.001)= 0; J21c(abs(J21c)<.001)= 0;J21d(abs(J21d)<.001)= 0;
J12a(abs(J12a)<.001)= 0; J12b(abs(J12b)<.001)= 0; J12c(abs(J12c)<.001)= 0;J12d(abs(J12d)<.001)= 0;
% calculate ground truth normals and metric tensor % x1 and x2
[Et,Ft,Gt,x1,x2,N,Pest]= create_ground_truth_normals_and_metric_tensor(Pgth,qgth,I1u,I1v,I2u,I2v);

tt = size(J21a,1);
% create E F G
E = (J21a.^2).*repmat(Et,tt,1) + 2.*J21a.*J21b.*repmat(Gt,tt,1) + (J21b.^2).*repmat(Ft,tt,1);
F = (J21c.^2).*repmat(Et,tt,1) + 2.*J21c.*J21d.*repmat(Gt,tt,1) + (J21d.^2).*repmat(Ft,tt,1);
G = (J21a.*J21c).*repmat(Et,tt,1) + (J21a.*J21d + J21b.*J21c).*repmat(Gt,tt,1) + (J21b.*J21d).*repmat(Ft,tt,1);



D1 = (J12a.*H21uua + J12c.*H21uub + J12b.*H21uva + J12d.*H21uvb)./3;
D2 = (J12a.*H21uva + J12c.*H21uvb + J12b.*H21vva + J12d.*H21vvb)./3;

k1 = J21a.*repmat(x1,tt,1) + J21b.*repmat(x2,tt,1) - D1;
k2 = J21c.*repmat(x1,tt,1) + J21d.*repmat(x2,tt,1) - D2;

% solve the equations

for i = 1: size(I2u,1) % for each input image
    Nest{i} = [k1(i,:);k2(i,:);1-I2u(i,:).*k1(i,:)-I2v(i,:).*k2(i,:)];
    Nest{i} = Nest{i}./repmat(sqrt(sum(Nest{i}.^2)),3,1);
    err_s{i} = acosd(dot(Nest{i},N{i+1}));
    beta{i} = nthroot((k1(i,:).^2 + k2(i,:).^2+(1-I2u(i,:).*k1(i,:)-I2v(i,:).*k2(i,:)).^2)./((x1.^2 + x2.^2+(1-I1u(1,:).*x1-I1v(1,:).*x2).^2).*Pest{i}(3,:).^4.*(J21a(i,:).*J21d(i,:)-J21b(i,:).*J21c(i,:)).^2),4);
    [Nest{i},idx]  = selectNormals( [I2u(i,:);I2v(i,:)], Nest{i}, Nest{i}, th);
    u_all = I2u(i,:); v_all = I2v(i,:);
    % find indices with no solution
    Nest{i}(:,idx) = []; u_all(:,idx) = []; v_all(:,idx) = []; beta{i}(:,idx)=[];
    % Integrate normals to find depth
    P_grid=calculate_depth(Nest{i},u_all,v_all,1e0);
    %P_grid=[u_all;v_all;ones(1,length(u_all))]./repmat(beta{i},3,1);
    % compare with ground truth
    [P2{i},err_p{i}] = compare_with_Pgth(P_grid,u_all,v_all,qgth{i+1},Pgth(3*i+1:3*i+3,:));
    [N{i},err_n{i}] = compare_with_Ngth(P2{i},qgth{i+1},Ngth(3*i+1:3*i+3,:));
end
