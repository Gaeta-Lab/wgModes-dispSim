%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script uses the WGMODES vector finite difference modesolver package
% to compute the variation of effective index of a waveguide with
% wavelength. WGMODES is a MATLAB-based open-source modesolver package
% written by Prof. Thomas E. Murphy's group at the University of Maryland.

% Web page: https://photonics.umd.edu/software/wgmodes/
% Please cite the following paper:
% A. B. Fallahkhair, K. S. Li and T. E. Murphy,
% "Vector Finite Difference Modesolver for Anisotropic Dielectric Waveguides",
% J. Lightwave Technol. 26(11), 1423-1431, (2008)
% 
% Created by Sai Kanth Dacha
% Last updated: August 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
close all
clc

mydir=pwd;
if ismac
    idcs=strfind(mydir,'/');
    parentdir=mydir(1:idcs(end)-1);
    destination_folder=[parentdir,'/'];
    funcdir=[mydir,'/Functions and dependencies/'];
elseif ispc
    idcs=strfind(mydir,'\');
    parentdir=mydir(1:idcs(end)-1);
    destination_folder=[parentdir,'\'];
    funcdir=[mydir,'\Functions and dependencies\'];
end
addpath(funcdir)

%% Physical constants
fundamental_constants;
hbar=physconsts.hbar; % in Joule seconds
epsilon0=physconsts.eps0;
c=physconsts.c;
mu0=1/(c^2*epsilon0);

%%
% Material options:
% 'sin': silicon nitride
% 'sio2': silica
% 'hto': high temperature oxide
% 'si': silicon
% 'mgf2': magnesium flouride
% 'diamond': diamond

% NOTE: All length dimensions here are in nm
lowCladMat='sio2';
h1 = 3000;           % SiO2 cladding height (referring to this as "height" assuming horizontal orientation of longer dimension of core)
coreMat='sin';
h2 = 730;           % SiN core height
upperCladMat='hto';
h3 = 2000;           % HTO cladding height

w = 1500;            % waveguide full-width
side = 2000;         % space on side of waveguide

dx=10; % Spatial grid size, in nm
dy=dx;
N_calc=4;

bend_radius_nm=100*10^3; % Ring bend radius, in nm

h_params.h1=h1; h_params.h2=h2; h_params.h3=h3;
h_params.w=w; h_params.side=side;
h_params.dx=dx; h_params.dy=dy;
h_params.bend_radius_nm=bend_radius_nm;
h_params.lowCladMat=lowCladMat; h_params.coreMat=coreMat; h_params.upperCladMat=upperCladMat;


lambda0_arr=500:50:3000; % Change this for start and stop wavelengths and step size (in nm)
neff_TEmode=lambda0_arr*0;
all_neffs=zeros(length(lambda0_arr),N_calc);

tPar=tic;
parfor (il=1:length(lambda0_arr),4)
    lambda0=lambda0_arr(il);
    neff=modecalc_fullgrid_step_func(h_params,lambda0,N_calc);
    neff_TEmode(il)=neff(1);
    all_neffs(il,:)=neff;
    il
end
toc(tPar)

%%
save([destination_folder,num2str(w),'by',num2str(h2),'_Rb',num2str(bend_radius_nm/10^3),'um_neff_vs_lam.mat'],'lambda0_arr','neff_TEmode','bend_radius_nm','all_neffs');
