%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script uses the WGMODES vector finite difference modesolver package
% to compute the variation of effective index of a waveguide with
% wavelength. WGMODES is a MATLAB-based open-source modesolver package
% written by Prof. Thomas E. Murphy's group at the University of Maryland.

% Web page: https://photonics.umd.edu/software/wgmodes/
% If you use this code, PLEASE CITE the following paper:
% A. B. Fallahkhair, K. S. Li and T. E. Murphy,
% "Vector Finite Difference Modesolver for Anisotropic Dielectric Waveguides",
% J. Lightwave Technol. 26(11), 1423-1431, (2008)
% 
% Created by Sai Kanth Dacha (saikanth.dacha@gmail.com)
% Last updated: May 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPLANATION:
% This code computes dispersion of various orders (beta0-beta5) at a desired wavelength
% for the fundamental TE spatial mode a SiN-core waveguide with a desired waveguide width, height, and bending radius.
% The dispersion is computed by (spline) interpolating a pre-computed multi-dimensional matrix of 
% effective indices.
% 
% Spatial grid size: 10 nm
% Bottom cladding: SiO2 (3 um thick)
% Top cladding: HTO (2 um thick)
% Core: SiN
% Sellmeier for SiN:
% eV_to_nm = 1239.842;
% out = sqrt(1+253.86/(9.1609^2-(eV_to_nm/wl_nm)^2)+0.040314/(0.001^2-(eV_to_nm/wl_nm)^2));
% 
% The pre-computed matrix is generated using the modesolver package referenced above for
% the following parameters:
% Waveguide width: 800 nm to 3000 nm in steps of 100 nm
% Waveguide height: 600 nm to 800 nm in steps of 40 nm
% Waveguide bending radius: [20,40,60,100,150,200,350,500,800,1500,4000,10000] um
% Wavelength: 400 nm to 2200 nm in steps of 40 nm
% 
% DEPENDENCIES:
% The following files must be in the same directory as this function:
% fullSweep_TE0_neff_vs_lam.mat
% calc_betas.m
% finite_diff_deriv.m

% INPUT:
% wg_width_nm: waveguide width (in nm) [between 900 nm and 3000 nm]
% wg_height_nm: waveguide height (in nm) [between 600 nm and 800 nm]
% wg_Rb_um: waveguide bending radius (in um) [between 20 um and 10000 um]
% selWl_nm: wavelength of interest (in nm) [between 500 nm and 2600 nm]
% 
% OUTPUT:
% betas_interp: MATLAB struct variable, with beta0-beta5 as the six fields

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function betas_interp=betas_interpFunc_gaetaLab(wg_width_nm,wg_height_nm,wg_Rb_um,selWl_nm)

load('fullSweep_TE0_neff_vs_lam.mat');

%% Physical constants
c=299792458;
physconsts.c=c;

wg_Rb_nm=wg_Rb_um*1000;

%%
beta0_fullSweep=zeros(length(h2_sweep),length(w_sweep),length(bend_radius_nm_sweep));
beta1_fullSweep=zeros(size(beta0_fullSweep));
beta2_fullSweep=zeros(size(beta0_fullSweep));
beta3_fullSweep=zeros(size(beta0_fullSweep));
beta4_fullSweep=zeros(size(beta0_fullSweep));
beta5_fullSweep=zeros(size(beta0_fullSweep));

for ih2=1:length(h2_sweep)
    for iw=1:length(w_sweep)
        for ib=1:length(bend_radius_nm_sweep)
            neff_isweep=neffs_fullSweep(:,ih2,iw,ib);
            betas_sim=calc_betas(neff_isweep',lambda0_arr,selWl_nm,physconsts);
            beta0_fullSweep(ih2,iw,ib)=betas_sim.beta0_laminterp(betas_sim.selWl_nm_ind);
            beta1_fullSweep(ih2,iw,ib)=betas_sim.beta1_laminterp(betas_sim.selWl_nm_ind);
            beta2_fullSweep(ih2,iw,ib)=betas_sim.beta2_laminterp(betas_sim.selWl_nm_ind);
            beta3_fullSweep(ih2,iw,ib)=betas_sim.beta3_laminterp(betas_sim.selWl_nm_ind);
            beta4_fullSweep(ih2,iw,ib)=betas_sim.beta4_laminterp(betas_sim.selWl_nm_ind);
            beta5_fullSweep(ih2,iw,ib)=betas_sim.beta5_laminterp(betas_sim.selWl_nm_ind);
        end
    end
end


fprintf(['\n' ...
    '\n' ...
    'For accurate results, please input values within the following limits:\n' ...
    '600 nm < height < 800 nm\n' ...
    '800 nm < width < 3000 nm\n' ...
    '20 um < bend radius < 10 mm\n']);

%%
bendRadMax_idx=length(bend_radius_nm_sweep);
bendRadMax=bend_radius_nm_sweep(bendRadMax_idx);

if wg_Rb_nm>bendRadMax
    fprintf('WARNING: Chosen bend radius > simulated maximum;\nDisplaying dispersion for simulated maximum (10 mm) \n');
    beta0_interp=interpn(h2_sweep,w_sweep,beta0_fullSweep(:,:,bendRadMax_idx),wg_height_nm,wg_width_nm,'spline');
    beta1_interp=interpn(h2_sweep,w_sweep,beta1_fullSweep(:,:,bendRadMax_idx),wg_height_nm,wg_width_nm,'spline');
    beta2_interp=interpn(h2_sweep,w_sweep,beta2_fullSweep(:,:,bendRadMax_idx),wg_height_nm,wg_width_nm,'spline');
    beta3_interp=interpn(h2_sweep,w_sweep,beta3_fullSweep(:,:,bendRadMax_idx),wg_height_nm,wg_width_nm,'spline');
    beta4_interp=interpn(h2_sweep,w_sweep,beta4_fullSweep(:,:,bendRadMax_idx),wg_height_nm,wg_width_nm,'spline');
    beta5_interp=interpn(h2_sweep,w_sweep,beta5_fullSweep(:,:,bendRadMax_idx),wg_height_nm,wg_width_nm,'spline');

elseif wg_Rb_nm<bend_radius_nm_sweep(1)
    fprintf('WARNING: Chosen bend radius < simulated minimum;\nDisplaying dispersion for simulated minimum (20 um) \n');
    beta0_interp=interpn(h2_sweep,w_sweep,beta0_fullSweep(:,:,1),wg_height_nm,wg_width_nm,'spline');
    beta1_interp=interpn(h2_sweep,w_sweep,beta1_fullSweep(:,:,1),wg_height_nm,wg_width_nm,'spline');
    beta2_interp=interpn(h2_sweep,w_sweep,beta2_fullSweep(:,:,1),wg_height_nm,wg_width_nm,'spline');
    beta3_interp=interpn(h2_sweep,w_sweep,beta3_fullSweep(:,:,1),wg_height_nm,wg_width_nm,'spline');
    beta4_interp=interpn(h2_sweep,w_sweep,beta4_fullSweep(:,:,1),wg_height_nm,wg_width_nm,'spline');
    beta5_interp=interpn(h2_sweep,w_sweep,beta5_fullSweep(:,:,1),wg_height_nm,wg_width_nm,'spline');
else
    beta0_interp=interpn(h2_sweep,w_sweep,bend_radius_nm_sweep,beta0_fullSweep,wg_height_nm,wg_width_nm,wg_Rb_nm,'spline');
    beta1_interp=interpn(h2_sweep,w_sweep,bend_radius_nm_sweep,beta1_fullSweep,wg_height_nm,wg_width_nm,wg_Rb_nm,'spline');
    beta2_interp=interpn(h2_sweep,w_sweep,bend_radius_nm_sweep,beta2_fullSweep,wg_height_nm,wg_width_nm,wg_Rb_nm,'spline');
    beta3_interp=interpn(h2_sweep,w_sweep,bend_radius_nm_sweep,beta3_fullSweep,wg_height_nm,wg_width_nm,wg_Rb_nm,'spline');
    beta4_interp=interpn(h2_sweep,w_sweep,bend_radius_nm_sweep,beta4_fullSweep,wg_height_nm,wg_width_nm,wg_Rb_nm,'spline');
    beta5_interp=interpn(h2_sweep,w_sweep,bend_radius_nm_sweep,beta5_fullSweep,wg_height_nm,wg_width_nm,wg_Rb_nm,'spline');
end


%%
betas_interp.beta0=beta0_interp;
betas_interp.beta1=beta1_interp;
betas_interp.beta2=beta2_interp;
betas_interp.beta3=beta3_interp;
betas_interp.beta4=beta4_interp;
betas_interp.beta5=beta5_interp;


