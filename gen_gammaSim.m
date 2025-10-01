%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script uses the WGMODES vector finite difference modesolver package
% to compute the spatial mode field profiles of a waveguide at a specified
% wavelength. WGMODES is a MATLAB-based open-source modesolver package
% written by Prof. Thomas E. Murphy's group at the University of Maryland.

% Web page: https://photonics.umd.edu/software/wgmodes/
% Please cite the following paper:
% A. B. Fallahkhair, K. S. Li and T. E. Murphy,
% "Vector Finite Difference Modesolver for Anisotropic Dielectric Waveguides",
% J. Lightwave Technol. 26(11), 1423-1431, (2008)
%
% Created by Sai Kanth Dacha
% Last updated: June 2023
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

%% Edit parameters in this section as desired

N_calc=4; % Number of modes to solve for
plot_mode_ind=1; % Mode index to be plotted (all modes are saved regardless)


% lambda0=1550; % Wavelength, in nm
% w_arr = 900:50:3000;            % waveguide full-width

lambda0_arr=900:50:2200; % Wavelength, in nm
w = 1500;            % waveguide full-width

side=2400;
dx=10; % Spatial grid size, in nm

bend_radius_nm=100*10^3; % Ring bend radius, in nm

n2=0.24*10^-18; % At 1550 nm, in m^2/W
%% Gamma Calculation
gammas_w_arr=zeros(N_calc,length(lambda0_arr));

for iw=1:length(lambda0_arr)
    [modes_EHfields,neff,spat_arrays]=modecalc_wg_fullgrid(w,lambda0_arr(iw),dx,N_calc,bend_radius_nm,side,plot_mode_ind,0,0);
    for in=1:length(neff)
        gammas_w_arr(in,iw)=calc_gam(modes_EHfields{in},n2,lambda0_arr(iw),spat_arrays,physconsts);
    end
    iw
end

% save([destination_folder,'gamma-vs-lam_',num2str(spat_arrays.h_params.h2),'_Rb',num2str(bend_radius_nm/10^3),'um_',num2str(lambda0),'nm.mat'],'lambda0_arr','w','spat_arrays','bend_radius_nm','gammas_w_arr');
