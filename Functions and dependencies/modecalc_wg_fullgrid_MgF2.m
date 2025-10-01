function [modes_EHfields,neff,spat_arrays]=modecalc_wg_fullgrid_MgF2(w,lambda0,dx,N_calc,bend_radius_nm,side,plot_mode_ind,cshift,plot_true)

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
coreMat='mgf2';
h2 = 1000;           % SiN core height
upperCladMat='hto';
h3 = 2000;           % HTO cladding height

% side = 2400;         % space on side of waveguide
dy=dx;

h_params.h1=h1;
h_params.h2=h2;
h_params.h3=h3;
h_params.w=w;

%% Index Grid
% n1 = sellmeier(lowCladMat,lambda0*10^-9);          % SiO2 lower cladding
n1 = 1;
n2 = sellmeier(coreMat,lambda0*10^-9);          % SiN core
% n3 = sellmeier(upperCladMat,lambda0*10^-9);          % HTO upper cladding
n3 = 1;

% fprintf (1,'generating index mesh...\n');
[x,y,xc,yc,nx,ny,eps,edges] = ...
    waveguidemeshfull([n1,n2,n3],[h1,h2,h3],h2,w/2,side,dx,dy); 

w2=100; h22 = 800;

eps(10:end-10,10:(h1-h22)/dx)=n2^2;
eps((side+w/2-w2/2)/dx:(side+w/2+w2/2)/dx,1:h1/dx)=n2^2;


% 
% [x,y,xc,yc,nx,ny,eps2,edges] = waveguidemeshfull([n1,n2,n3],[h1,h22,h3],h22,w2/2,side,dx,dy);
% 
% cshift2=0.5*(h22+h2)/dx;
% eps2=circshift(eps2,cshift2,2);
% 
% 
% [x,y,xc,yc,nx,ny,eps3,edges] = waveguidemeshfull([n1,n2,n3],[h1,h23,h3],h23,w3/2,side,dx,dy);
% cshift3=(0.5*(h23+h2)+h22)/dx;
% eps3=circshift(eps3,cshift3,2);
% 
% eps=eps1+eps2+eps3;


% Now we stretch out the mesh at the boundaries:
[x,y,xc,yc,dx_str,dy_str] = stretchmesh(x,y,[100,100,100,100],[4,4,4,4]);
n_grid=sqrt(eps);
% 
% % "Simulating" bending by adding a slope to the core index profile
% core_reg_logical=(n_grid>0.9*n2);
% bend_pad_width=400; % Width of region in cladding on either side of core that I simulate bending for
% bend_pad_ind=round(bend_pad_width/dx);
% 
% [co,ro,~]=find(core_reg_logical==1);
% hslice=unique(ro);
% wslice=unique(co');
% wslice=wslice(1)-bend_pad_ind:wslice(end)+bend_pad_ind;

% n_grid_bend=n_grid;
% 
% bend_factor=(1+(xc(wslice)-xc(wslice(round(length(wslice)/2))))/bend_radius_nm);
% bend_factor=bend_factor';
% 
% bend_matrix=1+0*hslice*wslice;
% bend_matrix=bend_matrix.*bend_factor; bend_matrix=bend_matrix';
% 
% n_grid_bend(wslice,hslice)=n_grid(wslice,hslice).*bend_matrix;
% 
% n_grid=circshift(n_grid,cshift,1);
% n_grid_bend=circshift(n_grid_bend,cshift,1);

eps=n_grid.^2;
% eps_bend=n_grid_bend.^2;

%% Mode Calculation

guess=max(max(n_grid));

boundary='0000'; % See wgmodes.m for details

t_start=tic();
[Hx,Hy,neff]=wgmodes(lambda0,guess,N_calc,dx_str,dy_str,eps,boundary);

neff=unique(abs(neff),'stable'); % Since we're not including scattering loss, this is okay
keep_modes_idx=find(neff>min([n1,n2,n3])); % Only keeping modes that would be guided under conventional TIR conditions
neff=neff(keep_modes_idx);
N_modes=length(neff);

toc(t_start);

modes_EHfields=cell(N_modes,1);
% modes_Hfields=cell(N_modes,1);

for i=1:N_modes
    ii=keep_modes_idx(i);
    hx=Hx(:,:,ii);
    hy=Hy(:,:,ii);
    [hz,ex,ey,ez]=postprocess(lambda0,neff(ii),hx,hy,dx_str,dy_str,eps,boundary);

    Efields_i(:,:,1)=ex;
    Efields_i(:,:,2)=ey;
    Efields_i(:,:,3)=ez;

    hx=hx(1:end-1,1:end-1);
    hy=hy(1:end-1,1:end-1);
    hz=hz(1:end-1,1:end-1);
    Hfields_i(:,:,1)=hx;
    Hfields_i(:,:,2)=hy;
    Hfields_i(:,:,3)=hz;

    mode_fields.E=Efields_i;
    mode_fields.H=Hfields_i;

    % Normalization method: \integ{E \times conj(H)}dA=1
    normfac=sqrt(sum(sum(dotprod_2dbeams(mode_fields,conj_mode_struct(mode_fields),1)))*dx*dx);
    
    mode_fields.E=Efields_i/normfac;
    mode_fields.H=Hfields_i/normfac;
    
    modes_EHfields{i}=mode_fields;
end

clear ex ey ez hx hy hz Efields_i Ex Ey Ez %Hfields_i Hx Hy Hz

[X,Y]=meshgrid(xc,flip(yc));
[Phi_grid,R_grid]=cart2pol(X,Y);
%%
if plot_true
    chosenMode_Efields=modes_EHfields{plot_mode_ind}; chosenMode_Efields=chosenMode_Efields.E;
    plot_mode_Efields(X,Y,chosenMode_Efields,'jet',h_params,0);
end

%% Exporting
spat_arrays.X=X;
spat_arrays.Y=Y;
spat_arrays.R_grid=R_grid;
spat_arrays.Phi_grid=Phi_grid;
spat_arrays.xc=xc;
spat_arrays.dx=dx;
spat_arrays.eps=eps;
spat_arrays.n_grid=n_grid;
% spat_arrays.n_grid_bend=n_grid_bend;
spat_arrays.h_params=h_params;
