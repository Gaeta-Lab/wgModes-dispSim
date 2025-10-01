function g=calc_gam(mode_EHfields,n2,lambda0,spat_arrays,physconsts)

hbar=physconsts.hbar; % in Joule seconds
epsilon0=physconsts.eps0;
c=physconsts.c;
mu0=1/(c^2*epsilon0);

w0=2*pi*c/(lambda0*10^-9);

%% Spatial arrays
dx=spat_arrays.dx;

%% Modes
mode_Efields=mode_EHfields.E;
A_eff_num=(sum(sum(sum(abs(mode_Efields).^2,3)))*dx*dx)^2;
A_eff_denom=sum(sum(sum(abs(mode_Efields).^4,3)))*dx*dx;

A_eff=A_eff_num/A_eff_denom * 10^-9*10^-9; % Multiplying twice by 10^-9 to account for the fact that the modesolver works in nm
g=n2*w0/(c*A_eff);

end
