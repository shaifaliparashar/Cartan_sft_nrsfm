
%ISOMETRIC REAL DATA TEST
clear all;
close all;

addpath(genpath('BBS'));
addpath(genpath('gloptipoly3'));
addpath(genpath('SeDuMi_1_3'));
addpath(genpath('schwarps'));
addpath('code');


load Kinect_paper.mat

th = 30;
n = 20;
idx = round(linspace(1,191,n)); % select 20 images to evaluate results on
idd = idx(2:end)-1;

H21uua = H21uua(idd,:);
H21uub = H21uub(idd,:);
H21uva = H21uva(idd,:);
H21uvb = H21uvb(idd,:);
H21vva = H21vva(idd,:);
H21vvb = H21vvb(idd,:);
I2u = I2u(idd,:);
I2v = I2v(idd,:);
detr = J21a(idd,:).*J21d(idd,:)-J21c(idd,:).*J21c(idd,:);
J21a = J21a(idd,:);
J21b = J21b(idd,:);
J21c = J21c(idd,:);
J21d = J21d(idd,:);
J12a = J12d(idd,:)./detr;
J12b = -J12b(idd,:)./detr;
J12c = -J12c(idd,:)./detr;
J12d = J12a(idd,:)./detr;

iddd = 1:5:1503;
for i = 1:length(idx)
    Pg(3*(i-1)+1:3*(i-1)+3,:) = Pgth(3*(idx(i)-1)+1:3*(idx(i)-1)+3,iddd);
    Ng(3*(i-1)+1:3*(i-1)+3,:) = Ngth(3*(idx(i)-1)+1:3*(idx(i)-1)+3,iddd);
    qg{i} = qgth{idx(i)}(1:2,iddd);
end

% ISOMETRIC/ CONFORMAL NRSFM
disp('solving isometric/conformal nrsfm')
[P_i,N_i,err_p_i,err_n_i] = isometric_nrsfm(Pg,Ng,qg,I1u,I1v,I2u,I2v,J21a,J21b,J21c,J21d,J12a,J12b,J12c,J12d,H21uva,H21uvb,1e-2);

% SKEWLESS NRSFM
disp('solving skewless nrsfm')
[P_s,N_s,err_p_s,err_n_s] = skewless_nrsfm(Pg,Ng,qg,I1u,I1v,I2u,I2v,J21a,J21b,J21c,J21d,J12a,J12b,J12c,J12d,H21uua,H21uub,H21uva,H21uvb,H21vva,H21vvb);


% ISOMETRIC/CONFORMAL SFT
disp('solving isometric/conformal sft')
[P23,Ni3,err_p3,err_n,P2c3,Nc,err_pc3,err_nc] = solve_iso_con_sft_final(Pg,Ng,qg,I1u,I1v,I2u,I2v,J21a,J21b,J21c,J21d,H21uua,H21uub,H21vva,H21vvb,th);

% SKEWLESS SFT
disp('solving skewless sft')
[P_sk,N_sk,err_p_sk,err_n_sk] = solve_skewless_sft(Pg,Ng,qg,I1u,I1v,I2u,I2v,J21a,J21b,J21c,J21d,J12a,J12b,J12c,J12d,H21uua,H21uub,H21uva,H21uvb,H21vva,H21vvb,th);

% EQUIAREAL SFT
disp('solving equiareal sft')
[P_e,N_e,err_p_e,err_n_e] = solve_equiareal_sft(Pg,Ng,qg,I1u,I1v,I2u,I2v,J21a,J21b,J21c,J21d,J12a,J12b,J12c,J12d,H21uua,H21uub,H21uva,H21uvb,H21vva,H21vvb,th);

% SMOOTH SFT
disp('solving smooth sft')
[P_sm,N_sm,err_p_sm,err_n_sm] = solve_smooth_sft(Pg,Ng,qg,I1u,I1v,I2u,I2v,J21a,J21b,J21c,J21d,J12a,J12b,J12c,J12d,H21uva,H21uvb,th);




