function [P2,N,err_p,err_n] = solve_smooth_sft(Pgth,Ngth,qgth,I1u,I1v,I2u,I2v,J21a,J21b,J21c,J21d,J12a,J12b,J12c,J12d,H21uva,H21uvb,th)
% calculate ground truth normals and metric tensor % x1 and x2
[Et,Ft,Gt,x1,x2,N,Pest]= create_ground_truth_normals_and_metric_tensor(Pgth,qgth,I1u,I1v,I2u,I2v);

tt = size(J21a,1);
% create E F G
k1 = J21a.*repmat(x1,tt,1) + J21b.*repmat(x2,tt,1) - (J12b.*H21uva + J12d.*H21uvb);
k2 = J21c.*repmat(x1,tt,1) + J21d.*repmat(x2,tt,1) - (J12a.*H21uva + J12c.*H21uvb);

% solve the equations

for i = 1: size(I2u,1) % for each input image
    Nest{i} = [k1(i,:);k2(i,:);1-I2u(i,:).*k1(i,:)-I2v(i,:).*k2(i,:)];
    Nest{i} = Nest{i}./repmat(sqrt(sum(Nest{i}.^2)),3,1);
    err_s{i} = min(acosd(dot(Nest{i},N{i+1})),acosd(dot(Nest{i},-N{i+1})));
    mean(err_s{i})
end
% solve the equations

for i = 1: size(I2u,1) % for each input image
    Nest{i} = [k1(i,:);k2(i,:);1-I2u(i,:).*k1(i,:)-I2v(i,:).*k2(i,:)];
    Nest{i} = Nest{i}./repmat(sqrt(sum(Nest{i}.^2)),3,1);
    err_s{i} = acosd(dot(Nest{i},N{i+1}));
    [Nest{i},idx]  = selectNormals( [I2u(i,:);I2v(i,:)], Nest{i}, Nest{i}, th);
    u_all = I2u(i,:); v_all = I2v(i,:);
    % find indices with no solution
    Nest{i}(:,idx) = []; u_all(:,idx) = []; v_all(:,idx) = []; 
    % Integrate normals to find depth
    P_grid=calculate_depth(Nest{i},u_all,v_all,1e-2);
    % compare with ground truth
    [P2{i},err_p{i}] = compare_with_Pgth(P_grid,u_all,v_all,qgth{i+1},Pgth(3*i+1:3*i+3,:));
    [N{i},err_n{i}] = compare_with_Ngth(P2{i},qgth{i+1},Ngth(3*i+1:3*i+3,:));
end
