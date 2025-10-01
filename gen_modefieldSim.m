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
lambda0=1550; % Wavelength, in nm

w = 2100;            % waveguide full-width
side=2400;
dx=10; % Spatial grid size, in nm

bend_radius_nm=inf; % Ring bend radius, in nm

%% Mode Calculation

[modes_EHfields,neff,spat_arrays]=modecalc_wg_fullgrid(w,lambda0,dx,N_calc,bend_radius_nm,side,plot_mode_ind,0,1);

% save([destination_folder,num2str(w),'by',num2str(h2),'_Rb',num2str(bend_radius_nm/10^3),'um_modeEfields_',num2str(lambda0),'nm.mat'],'lambda0','neff','spat_arrays','h_params','bend_radius_nm','modes_Efields');
