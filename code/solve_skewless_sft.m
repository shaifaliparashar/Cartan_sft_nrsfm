function [P2g,N,err_p,err_n] = solve_skewless_sft(Pgth,Ngth,qgth,I1u,I1v,I2u,I2v,J21a,J21b,J21c,J21d,J12a,J12b,J12c,J12d,H21uua,H21uub,H21uva,H21uvb,H21vva,H21vvb,th)
J21a(abs(J21a)<.001)= 0; J21b(abs(J21b)<.001)= 0; J21c(abs(J21c)<.001)= 0;J21d(abs(J21d)<.001)= 0;
J12a(abs(J12a)<.001)= 0; J12b(abs(J12b)<.001)= 0; J12c(abs(J12c)<.001)= 0;J12d(abs(J12d)<.001)= 0;
% calculate ground truth normals and metric tensor % x1 and x2
[Et,Ft,Gt,x1,x2,N,Pest]= create_ground_truth_normals_and_metric_tensor(Pgth,qgth,I1u,I1v,I2u,I2v);

tt = size(J21a,1);
% create E F G
E = (J21a.^2).*repmat(Et,tt,1) + 2.*J21a.*J21b.*repmat(Gt,tt,1) + (J21b.^2).*repmat(Ft,tt,1);
F = (J21c.^2).*repmat(Et,tt,1) + 2.*J21c.*J21d.*repmat(Gt,tt,1) + (J21d.^2).*repmat(Ft,tt,1);
G = (J21a.*J21c).*repmat(Et,tt,1) + (J21a.*J21d + J21b.*J21c).*repmat(Gt,tt,1) + (J21b.*J21d).*repmat(Ft,tt,1);

A1 = J21a.*repmat(Et,tt,1) + J21b.*repmat(Gt,tt,1);
A2 = J21a.*repmat(Gt,tt,1) + J21b.*repmat(Ft,tt,1);

C1 = J21c.*repmat(Et,tt,1) + J21d.*repmat(Gt,tt,1);
C2 = J21c.*repmat(Gt,tt,1) + J21d.*repmat(Ft,tt,1);

B1 = (C1.*H21uua + C2.*H21uub) + (A1.*H21uva + A2.*H21uvb);
B2 = (C1.*H21uva + C2.*H21uvb) + (A1.*H21vva + A2.*H21vvb);

D1 = (-J12a.*H21uua - J12c.*H21uub - J12b.*H21uva - J12d.*H21uvb)./(J21a.*J21d-J21b.*J21c);
D2 = (-J12a.*H21uva - J12c.*H21uvb - J12b.*H21vva - J12d.*H21vvb)./(J21a.*J21d-J21b.*J21c);



% eq1 = E*F*g12n*x2b + (B1 -E*G*(c*x1+d*x2) + (D1*F*E))*g22b = 0
P1 = (F.*D1 -E.*(J21c.*repmat(x1,tt,1)+J21d.*repmat(x2,tt,1))) + B1;
P2 = (F.*D2 -F.*(J21a.*repmat(x1,tt,1)+J21b.*repmat(x2,tt,1))) + B2;

% eq2 = E*F*g12n*x1b + (B2 -F*G*(a*x1+b*x2) + (D2*F*E))*g11b = 0
e = 1 + I2u.^2 + I2v.^2;


for i = 1: size(I2u,1) % for each input image
    
    for j = 1:size(I2v,2) % at each point
        
        %syms x1b x2b real
        mpol x1b x2b real
        E1 = [G(i,j)*e(i,j),P1(i,j)*e(i,j)-2*I2u(i,j)*G(i,j),-P1(i,j)*I2v(i,j),G(i,j)-I2u(i,j)*P1(i,j)];
        E2 = [G(i,j)*e(i,j),P2(i,j)*e(i,j)-2*I2v(i,j)*G(i,j),-P2(i,j)*I2u(i,j),G(i,j)-I2u(i,j)*P2(i,j)];
        %S = solve(E1*[x1b*x2b^2;x2b*x1b;x2b^2;x2b] == -P2(i,j)-P4(i,j), E1*[x1b*x2b^2;x2b*x1b;x2b^2;x2b] == -P3(i,j)-P4(i,j),'real',true);
        pol = (E1*[x2b*x1b^2;x2b*x1b;x1b;x2b])^2 + (E2*[x2b^2*x1b;x2b*x1b;x2b;x1b])^2;
        P = msdp(min(pol));
        [status,obj,M] = msol(P);
        if status ==1
            res = double([x1b x2b]);
        else
            res = [0,0];
        end
        k11(i,j) = res(1);
        k21(i,j) = res(2);

    end
end


for i = 1: size(I2u,1) % for each input image
    Nest1{i} = [k11(i,:);k21(i,:);1-I2u(i,:).*k11(i,:)-I2v(i,:).*k21(i,:)];
    Nest1{i} = Nest1{i}./repmat(sqrt(sum(Nest1{i}.^2)),3,1);
    [Nest{i},idx]  = selectNormals( [I2u(i,:);I2v(i,:)], Nest1{i}, Nest1{i}, th);
    u_all = I2u(i,:); v_all = I2v(i,:);
    % find indices with no solution
    Nest{i}(:,idx) = []; u_all(:,idx) = []; v_all(:,idx) = [];
    % Integrate normals to find depth
    P_grid=calculate_depth(Nest{i},u_all,v_all,1e0);
    % compare with ground truth
    [P2g{i},err_p{i}] = compare_with_Pgth(P_grid,u_all,v_all,qgth{i+1},Pgth(3*i+1:3*i+3,:));
    [N{i},err_n{i}] = compare_with_Ngth(P2g{i},qgth{i+1},Ngth(3*i+1:3*i+3,:));
    
end
