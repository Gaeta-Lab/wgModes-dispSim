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
addpath(destination_folder)
addpath(funcdir)

filename_str='2100by690_Rb100um_neff_vs_lam';
load([filename_str,'.mat']);
selModeInd=1;

%%
w_ind_end=strfind(filename_str,'by')-1;
wg_w=str2double(filename_str(1:w_ind_end));
h_ind_st=strfind(filename_str,'by')+2;
wg_h=str2double(filename_str(h_ind_st:h_ind_st+2));
Rb=str2double(filename_str(h_ind_st+6:h_ind_st+8));
% Rb=str2double(filename_str(h_ind_st+6:h_ind_st+8));

%% Physical constants
fundamental_constants;
hbar=physconsts.hbar; % in Joule seconds
epsilon0=physconsts.eps0;
c=physconsts.c;
mu0=1/(c^2*epsilon0);

%% Calculating derivative of neff with respect to wavelength

neff_selMode=neff_TEmode;
% neff_selMode=all_neffs(:,selModeInd)';

lambda_arr=lambda0_arr*10^-9;
w_arr=2*pi*c./lambda_arr;

lambda_interp=(600:1:2100)*10^-9;
w_interp=min(2*pi*c./lambda_interp):2*pi*100*10^9:max(2*pi*c./lambda_interp);

beta0_calc=w_arr.*neff_selMode/c;
beta0_laminterp=interp1(lambda_arr,beta0_calc,lambda_interp);
beta0_omeginterp=interp1(w_arr,beta0_calc,w_interp);

dlambda=lambda_arr(2)-lambda_arr(1);

dneff_dlambda=finite_diff_deriv(neff_selMode,dlambda); dneff_dlambda=dneff_dlambda';
d2neff_dlamba2=finite_diff_deriv(dneff_dlambda,dlambda); d2neff_dlamba2=d2neff_dlamba2';


%% Calculating beta1 and interpolating

n_g=neff_selMode-lambda_arr.*dneff_dlambda;
beta1_calc=n_g/c;
L_rt=700.461*10^-6;
fsr=1./(beta1_calc*L_rt);

n_eff_laminterp=interp1(lambda_arr,neff_selMode,lambda_interp);
n_eff_omeginterp=interp1(w_arr,neff_selMode,w_interp);

n_g_laminterp=interp1(lambda_arr,n_g,lambda_interp);
n_g_omeginterp=interp1(w_arr,n_g,w_interp);

beta1_laminterp=interp1(lambda_arr,beta1_calc,lambda_interp);
beta1_omeginterp=interp1(w_arr,beta1_calc,w_interp);

fsr_interp=interp1(lambda_arr,fsr,lambda_interp);
[~,lambda_1550_ind]=min(abs(lambda_interp-1550*10^-9));


%% Calculating beta2 and interpolating
dbeta1_dlambda=finite_diff_deriv(beta1_calc,dlambda)';
beta2_calc=-lambda_arr.^2.*dbeta1_dlambda/(2*pi*c);
% D=-lambda0_arr.*d2neff_dlamba2/c;

beta2_laminterp=interp1(lambda_arr,beta2_calc,lambda_interp);
beta2_omeginterp=interp1(w_arr,beta2_calc,w_interp);


%% Calculating beta3 and interpolating
dbeta2_dlambda=finite_diff_deriv(beta2_calc,dlambda)';
beta3_calc=-lambda_arr.^2.*dbeta2_dlambda/(2*pi*c);

beta3_laminterp=interp1(lambda_arr,beta3_calc,lambda_interp);
beta3_omeginterp=interp1(w_arr,beta3_calc,w_interp);

%% Calculating beta4 and interpolating
dbeta3_dlambda=finite_diff_deriv(beta3_calc,dlambda)';
beta4_calc=-lambda_arr.^2.*dbeta3_dlambda/(2*pi*c);

beta4_laminterp=interp1(lambda_arr,beta4_calc,lambda_interp);
beta4_omeginterp=interp1(w_arr,beta4_calc,w_interp);

%% Calculating beta5 and interpolating
dbeta4_dlambda=finite_diff_deriv(beta4_calc,dlambda)';
beta5_calc=-lambda_arr.^2.*dbeta4_dlambda/(2*pi*c);

beta5_laminterp=interp1(lambda_arr,beta5_calc,lambda_interp);
beta5_omeginterp=interp1(w_arr,beta5_calc,w_interp);

%% Plotting

figure(1)
plot(lambda_arr*10^9,beta2_calc*10^27)
xlabel('Wavelength (nm)'); ylabel('Calculated \beta_2 (ps^2/km)')
yyaxis right; plot(lambda_arr*10^9,fsr/10^9); hold on;
plot(lambda_arr*10^9,lambda_arr*0+fsr_interp(lambda_1550_ind)/10^9,'--')
ylabel('FSR (GHz)')
set(gca,'fontsize',16); grid on;

title(['Simulated GVD curve for ',num2str(wg_w),' nm x ',num2str(wg_h),' nm SiN/SiO2 bent waveguide (R_b = ',num2str(Rb),' \mum)'])
set(gcf,'Position',[400, 300, 1000, 500]);

fprintf(['beta_2 at 1550 nm = ',num2str(round(10^27*beta2_laminterp(lambda_1550_ind),2)),' ps2/km','\n'])



%% Saving

betas_sim.beta0_calc=beta0_calc;
betas_sim.beta1_calc=beta1_calc;
betas_sim.beta2_calc=beta2_calc;
betas_sim.beta3_calc=beta3_calc;
betas_sim.beta4_calc=beta4_calc;
betas_sim.beta5_calc=beta5_calc;
betas_sim.lambda_arr=lambda_arr;

betas_sim.beta0_laminterp=beta0_laminterp;
betas_sim.beta1_laminterp=beta1_laminterp;
betas_sim.beta2_laminterp=beta2_laminterp;
betas_sim.beta3_laminterp=beta3_laminterp;
betas_sim.beta4_laminterp=beta4_laminterp;
betas_sim.beta5_laminterp=beta5_laminterp;
betas_sim.lambda_interp=lambda_interp;

betas_sim.beta0_omeginterp=beta0_omeginterp;
betas_sim.beta1_omeginterp=beta1_omeginterp;
betas_sim.beta2_omeginterp=beta2_omeginterp;
betas_sim.beta3_omeginterp=beta3_omeginterp;
betas_sim.beta4_omeginterp=beta4_omeginterp;
betas_sim.beta5_omeginterp=beta5_omeginterp;
betas_sim.w_interp=w_interp;

betas_sim.n_g=n_g;
betas_sim.n_g_laminterp=n_g_laminterp;
betas_sim.n_g_omeginterp=n_g_omeginterp;

betas_sim.neff=neff_selMode;
betas_sim.neff_laminterp=n_eff_laminterp;
betas_sim.neff_omeginterp=n_eff_omeginterp;

% betas_sim.n_grid=n_grid;
% betas_sim.n_grid_bend=n_grid_bend;
betas_sim.bend_radius_nm=bend_radius_nm;

% save([destination_folder,num2str(wg_w),'by',num2str(wg_h),'_Rb',num2str(Rb),'um_betas_vs_lam.mat'],'betas_sim');