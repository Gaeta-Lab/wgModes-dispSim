function betas_sim=calc_betas(neff_selMode,lambda0_arr,selWl_nm,physconsts)

c=physconsts.c;

%% Calculating derivative of neff with respect to wavelength

lambda_arr=lambda0_arr*10^-9;
w_arr=2*pi*c./lambda_arr;

lambda_interp=(lambda0_arr(1):1:lambda0_arr(end))*10^-9;
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

n_eff_laminterp=interp1(lambda_arr,neff_selMode,lambda_interp);
n_eff_omeginterp=interp1(w_arr,neff_selMode,w_interp);

n_g_laminterp=interp1(lambda_arr,n_g,lambda_interp);
n_g_omeginterp=interp1(w_arr,n_g,w_interp);

beta1_laminterp=interp1(lambda_arr,beta1_calc,lambda_interp);
beta1_omeginterp=interp1(w_arr,beta1_calc,w_interp);

[~,selWl_nm_ind]=min(abs(lambda_interp-selWl_nm*10^-9));


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

%%
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

betas_sim.selWl_nm_ind=selWl_nm_ind;