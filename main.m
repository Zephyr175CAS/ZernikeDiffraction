clear; 
clc; 
close all; 
rng(2025);

lambda = 532e-9; 
f = 1; 
z = 0.5; 
R = 2e-3;
N = 512; 
overs = 2;

M=12; 
T=8; 
waves_sigma=0.25;
coeffs=(2*pi)*(waves_sigma*randn(M,T)); 
coeffs(:,1)=0;

out = simulate_zernike_imaging('f',f,'z',z,'R',R,'lambda',lambda,'N',N,'oversample',overs,'coeffs',coeffs);

% 3x4布局
fig1 = plot_wavefront_grid(out, struct('title',"入瞳相位","wrap",true,'layout',[3,4],'fontSize',12));
fig2 = plot_intensity_grid(out, struct('title',"z处强度(归一化)",'layout',[3,4],'fontSize',12));

% exportgraphics(fig1,'wavefront_grid.pdf','ContentType','vector');
% exportgraphics(fig2,'intensity_grid.pdf','ContentType','vector');
