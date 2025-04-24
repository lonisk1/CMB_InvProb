%% Driver_Ex1.m
%
% This script compares the reonstruction methods for the iterative
% reconstruction methods for imaging cosmic strings using
% Cosmic Microwave Background (CMB) data. This simulates Example #1.
%
% Chung, Onisk, and Wang. "ITERATIVE RECONSTRUCTION METHODS FOR
%   COSMOLOGICAL X-RAY TOMOGRAPHY", 2025
%
% Note: This is a demonstration of the code and does not match the results 
% presented in the paper.

clear
clc
close all

%% Sparse matrix for low resolution grid
% Detector Grid Definition
T_model = 4; % height of detectors. This value comes from model, f, in paper

% detector grid is the same xy as the image
n = 51; % this is an arbitrary number for grid-ing
range = 3+T_model; % Comes from spatial bounds from model of string (for us this is [-3,3])
xvec = linspace(-range,range,n);
yvec = linspace(-range,range,n);

% Scale up to integer grid
scale_d = 2*range/(n-1); %scaling param determined s.t [-range,range] with n points is scaled to
%                         [-n/2,n/2] integer grid with n points
xvec = round(xvec/scale_d);
yvec = round(yvec/scale_d);

[x_d,y_d] = ndgrid(xvec, yvec); % detector grid (integer) centered at zero
n_d = size(x_d,1);

% Source Grid Definition
source_rng = range; % radius of cone plus detector limit
scaled_rng_s = round(source_rng/scale_d);

xvec_s = -scaled_rng_s : scaled_rng_s;
yvec_s = -scaled_rng_s : scaled_rng_s;
[x_s,y_s] = ndgrid(xvec_s, yvec_s);

% Temporal Grid Definition
temporalSlices = 20; %arbitrarily chosen
rad = T_model; % non-scaled radius from model: should be T_model since require cone triangle to be 45-45-90
offset = (T_model/temporalSlices)/2; % used so that lowest/highest plane don't sit on source/detector grid
Emin = offset; Emax = T_model-offset; n_E = temporalSlices; % set height of detector and temporal slices
E = linspace(Emin,Emax,n_E); % non-scaled temporal discretization
E = E/scale_d; % scaled temporal discret. to match integer grid
T = T_model/scale_d; % scaled distance between source and det.: to match integer grid
n_g = (n-1)/2;
rad = rad/scale_d; % scaled radius (this info is duplicated from T above)

visual = 0; % only save the data from the detectors
tic,
A_lr = getRayTraceMat(x_s, y_s, x_d, y_d, rad, T, E, n, n_E, n_g,visual);
toc

%% Get true string
params.example = 1; % Example #1 = string in motion;
params.a = 0.05;
params.c = 1.5; %either 0,1,1.5 for Ex.1 in paper

scale = 2; % high-resolution to low-resolution scale
n_E = 20; % number of time points for low-resolution scale
xtrue = getstring(scale*n-1, n_E, params);

% Get LR string
[X,Y] = meshgrid(-(n-1):(n-1),-(n-1):(n-1));
[Xi,Yi] = meshgrid(-(n-1):scale:(n-1),-(n-1):scale:(n-1));
xtrue_lr = zeros(n,n,n_E);
for i = 1:size(xtrue,3)
  xtrue_lr(:,:,i) = interp2(X,Y,xtrue(:,:,i),Xi,Yi);
end

%% Sparse matrix for high resolution grid
% Temporal Grid Definition
temporalSlices = scale*n_E; % want DOUBLE the low-res number!
rad = T_model; % non-scaled radius from model: should be T_model since require cone triangle to be 45-45-90
offset = (T_model/temporalSlices)/2; % used so that lowest/highest plane don't sit on source/detector grid
Emin = offset; Emax = T_model-offset; n_E = temporalSlices; % set height of detector and temporal slices
E = linspace(Emin,Emax,n_E); % non-scaled temporal discretization
E = E/scale_d; % scaled temporal discret. to match integer grid
T = T_model/scale_d; % scaled distance between source and det.: to match integer grid
n_g = (n-1)/2;
rad = rad/scale_d; % scaled radius (this info is duplicated from T above)
tic
A_hr = getRayTraceMat2(x_s, y_s, x_d, y_d, rad, T, E, n, n_E, n_g, visual, scale);
toc

%% Generate observations using the HR grid
xtrue = getstring(scale*n-1, n_E, params);
b = (1/scale)*(A_hr*xtrue(:)); % here we scale b_hr down since we 2x number of time slices

%% Set-up Numerical Experiments (add noise)
n_E = 20; % want to ensure that we get the LR grid
xtrue = xtrue_lr; % want the coarse truth here
A = A_lr; % want to use coarse A for reconstructions
nlevel = 0.03; % noise level
etaLevel = 1.5; % parameter for Discrepancy Principle
bn = b + WhiteNoise(b,nlevel,7); % 7 is random number seed

opt.x_true = xtrue(:);
opt.MaxIter = 400;
opt.RegParam = 'discrep';
opt.eta = 1.01;
opt.NoiseLevel = nlevel;

%% Landweber Reconstruction
[X_land,info_land,ext_info] = landweber(A,bn,1:opt.MaxIter,zeros(size(xtrue(:))));
norm_b = norm(bn);
for i = 1:size(X_land,2)
  err_land(i) = norm(X_land(:,i) - xtrue(:)) / norm(xtrue(:));
  rsnrm_land(i) = norm(A*X_land(:,i) - bn) / norm_b;
end

% Determine what iter. DP satisfied for Landweber
breakout = etaLevel*nlevel;
for i = 1:length(rsnrm_land)
  if rsnrm_land(i) <= breakout
    landBreak_idx = i;
    break
  end
  if i == length(rsnrm_land)
    warning('Landweber did not reach breakout criterion')
    landBreak_idx = 400;
  end
end
[min_RRE_Land,idx_min_Land] = min(err_land); %determine idx of min error
x_land = reshape(X_land(:,idx_min_Land),n,n,n_E); %best RRE

%% Tikhonov Reconstruction
solver = 'tikhonov';
opt.NoStop = 'on';
opt.x_true = xtrue(:);
opt.MaxIter = 400;
[x_Tik, info_Tik] = IRhybrid_lsqr(A, bn(:), 1:opt.MaxIter,opt);
rsnrm_Tik = info_Tik.Rnrm;

breakout = etaLevel*nlevel;
for i = 1:length(rsnrm_Tik)
  if rsnrm_Tik(i) <= breakout
    TikBreak_idx = i;
    break
  end
  if i == length(rsnrm_Tik)
    warning('Tik did not reach breakout criterion')
    TikBreak_idx = 400;
  end
end
[min_RRE_Tik,idx_min_Tik] = min(info_Tik.Enrm);
x_Tik = reshape(x_Tik(:,idx_min_Tik),n,n,n_E);

%% FISTA Reconstruction
opt.MaxIter = 400;
opt.NoStop = 'on';

lambda = 1.5875; % Ex#1_c1.5_a0.05 (updated 3/1/25)
opt.RegParam = lambda;
[x_fista,info_fista] = IRfista(A,bn,opt);
x_fista = reshape(x_fista,n,n,n_E);

breakout = etaLevel*nlevel;
for i = 1:length(info_fista.Enrm)
  if info_fista.Rnrm(i) <= breakout
    fistaBreak_idx = i;
    break
  end
  if i == length(info_fista.Enrm)
    warning('fista did not reach breakout criterion')
    fistaBreak_idx = 400;
  end
end
[min_RRE_fista,idx_min_fista] = min(info_fista.Enrm);

%% Relative errors
figure(8), hold off
plot(err_land,'--b','Linewidth',1.4), hold on
plot(info_Tik.Enrm,'-r','Linewidth',1.4)
plot(info_fista.Enrm,'-.m','Linewidth',1.4)
plot(idx_min_Land,err_land(idx_min_Land),'bo','Markersize',8,'Linewidth',1.4) % for Land min Error
plot(landBreak_idx,err_land(landBreak_idx),'b*','Markersize',8,'Linewidth',1.4) % for Land breakout Error

plot(idx_min_fista,info_fista.Enrm(idx_min_fista),'mo','Markersize',8,'Linewidth',1.4) % for FISTA min Error
plot(fistaBreak_idx,info_fista.Enrm(fistaBreak_idx),'m*','Markersize',8,'Linewidth',1.4) % for FISTA breakout Error

plot(idx_min_Tik,info_Tik.Enrm(idx_min_Tik),'ro','Markersize',8,'Linewidth',1.4) % for Tik min Error
plot(TikBreak_idx,info_Tik.Enrm(TikBreak_idx),'r*','Markersize',8,'Linewidth',1.4) % for Tik breakout Error

set(gca,'fontsize',14)
hold off
lgd = legend('Landweber','Tikhonov','FISTA')
fontsize(lgd,20,'points')
axis([0 max([length(info_Tik.Enrm) length(err_land)]) 0 0.45])
xlabel('Iter. No.')
ylabel('RRE')

%% Relative Residuals
breakLvl = etaLevel*nlevel;
noiseVec = breakLvl.*ones(400,1);
figure(9)
semilogy(rsnrm_land,'--b','Linewidth',1.4), hold on
semilogy(rsnrm_Tik,'-r','Linewidth',1.4)
semilogy(info_fista.Rnrm,'-.m','Linewidth',1.4)
semilogy(noiseVec,':k','Linewidth',1.4)
set(gca,'fontsize',14)
hold off
lgd = legend('Landweber','Tikhonov','FISTA')
fontsize(lgd,20,'points')
axis([0 max([length(info_Tik.Rnrm) length(rsnrm_land)]) 0.29e-1 0.5e-1])
xlabel('Iter. No.')
ylabel('RRN')

%% All images for same time point

figure(10),
name = gray;
t_idx = 2;
subplot(3,4,1), imshow(xtrue(:,:,t_idx),[],'border','tight',colormap=name)
ylabel(sprintf('t=%d',t_idx), 'FontSize', 15,'Rotation',0,'VerticalAlignment','middle');
subplot(3,4,2), imshow(x_fista(:,:,t_idx),[],'border','tight',colormap=name),
subplot(3,4,3), imshow(x_Tik(:,:,t_idx),[],'border','tight',colormap=name),
subplot(3,4,4), imshow(x_land(:,:,t_idx),[],'border','tight',colormap=name),

t_idx = 10;
subplot(3,4,5), imshow(xtrue(:,:,t_idx),[],'border','tight',colormap=name),
ylabel(sprintf('t=%d',t_idx), 'FontSize', 15,'Rotation',0,'VerticalAlignment','middle');
subplot(3,4,6), imshow(x_fista(:,:,t_idx),[],'border','tight',colormap=name),
subplot(3,4,7), imshow(x_Tik(:,:,t_idx),[],'border','tight',colormap=name),
subplot(3,4,8), imshow(x_land(:,:,t_idx),[],'border','tight',colormap=name),

t_idx = 20;
subplot(3,4,9), imshow(xtrue(:,:,t_idx),[],'border','tight',colormap=name),
ylabel(sprintf('t=%d',t_idx), 'FontSize', 15,'Rotation',0,'VerticalAlignment','middle');
axes_handle = gca;
axes_handle.XLabel.Position = [26.0000 55 1.0000];
xlabel('True','FontSize', 14)
subplot(3,4,10), imshow(x_fista(:,:,t_idx),[],'border','tight',colormap=name),
axes_handle = gca;
axes_handle.XLabel.Position(2) = 55;
xlabel('FISTA','FontSize', 14)
subplot(3,4,11), imshow(x_Tik(:,:,t_idx),[],'border','tight',colormap=name),
axes_handle = gca;
axes_handle.XLabel.Position(2) = 55;
xlabel('Tikhonov','FontSize', 14)
subplot(3,4,12), imshow(x_land(:,:,t_idx),[],'border','tight',colormap=name),
axes_handle = gca;
axes_handle.XLabel.Position(2) = 55;
xlabel('Landweber','FontSize', 14)