function [modes_EHfields,neff,spat_arrays]=modecalc_ringbus_fullgrid(w_bus,w_ring,gap,lambda0,dx,N_calc,bend_radius_bus_nm,bend_radius_ring_nm,side,plot_mode_ind,plot_true)

% NOTE: All length dimensions here are in nm
lowCladMat='sio2';
h1 = 3000;           % SiO2 cladding height (referring to this as "height" assuming horizontal orientation of longer dimension of core)
coreMat='sin';
h2 = 730;           % SiN core height
upperCladMat='hto';
h3 = 2000;           % HTO cladding height

% side = 5000;         % space on side of waveguides
dy=dx;

h_params.h1=h1;
h_params.h2=h2;
h_params.h3=h3;
h_params.w_r=w_ring;
h_params.w_b=w_bus;
h_params.gap=gap;

%% Generating index grid (Bus)

n1 = sellmeier(lowCladMat,lambda0*10^-9);
n2 = sellmeier(coreMat,lambda0*10^-9);
n3 = sellmeier(upperCladMat,lambda0*10^-9);

fprintf (1,'generating index mesh...\n');
[x,y,xc,yc,nx,ny,eps_bus,edges_bus] = ...
    waveguidemeshfull([n1,n2,n3],[h1,h2,h3],h2,w_bus/2,side,dx,dy); 
n_grid_bus=sqrt(eps_bus);

% "Simulating" bending by adding a slope to the core index profile
core_reg_logical_bus=(n_grid_bus>0.9*n2);
bend_pad_width=400; % Width of region in cladding on either side of core that I simulate bending for
bend_pad_ind=round(bend_pad_width/dx);

[co,ro,~]=find(core_reg_logical_bus==1);
hsliceBus=unique(ro);
wsliceBus=unique(co');
wsliceBus=wsliceBus(1)-bend_pad_ind:wsliceBus(end)+bend_pad_ind;

core_reg_logical_bus(wsliceBus,hsliceBus)=1;

n_grid_bend_bus=n_grid_bus;

bend_factor=(1+(xc(wsliceBus)-xc(wsliceBus(round(length(wsliceBus)/2))))/bend_radius_bus_nm);
bend_factor=bend_factor';

bend_matrix=1+0*hsliceBus*wsliceBus;
bend_matrix=bend_matrix.*bend_factor; bend_matrix=bend_matrix';

n_grid_bend_bus(wsliceBus,hsliceBus)=n_grid_bus(wsliceBus,hsliceBus).*bend_matrix;
eps_bend_bus=n_grid_bend_bus.^2;

%% Generating index grid (Ring)

fprintf (1,'generating index mesh...\n');
[x,y,xc,yc,nx,ny,eps_ring,edges_ring] = ...
    waveguidemeshfull([n1,n2,n3],[h1,h2,h3],h2,w_ring/2,side,dx,dy); 
n_grid_ring=sqrt(eps_ring);


% "Simulating" bending by adding a slope to the core index profile
core_reg_logical_ring=(n_grid_ring>0.9*n2);
bend_pad_width=400; % Width of region in cladding on either side of core that I simulate bending for

if gap<=bend_pad_width
    bend_pad_width=gap/2;
end

bend_pad_ind=round(bend_pad_width/dx);

[co,ro,~]=find(core_reg_logical_ring==1);
hsliceRing=unique(ro);
wsliceRing=unique(co');
wsliceRing=wsliceRing(1)-bend_pad_ind:wsliceRing(end)+bend_pad_ind;

core_reg_logical_ring(wsliceRing,hsliceRing)=1;

n_grid_bend_ring=n_grid_ring;

bend_factor=(1+(xc(wsliceRing)-xc(wsliceRing(round(length(wsliceRing)/2))))/bend_radius_ring_nm);
bend_factor=bend_factor';

bend_matrix=1+0*hsliceRing*wsliceRing;
bend_matrix=bend_matrix.*bend_factor; bend_matrix=bend_matrix';

n_grid_bend_ring(wsliceRing,hsliceRing)=n_grid_ring(wsliceRing,hsliceRing).*bend_matrix;
eps_bend_ring=n_grid_bend_ring.^2;

%% Combining bus and ring index grids
n_grid_ring=circshift(n_grid_ring,-(gap+0.5*(w_ring+w_bus))/dx,1);
n_grid_bend_ring=circshift(n_grid_bend_ring,-(gap+0.5*(w_ring+w_bus))/dx,1);
core_reg_logical_ring=circshift(core_reg_logical_ring,-(gap+0.5*(w_ring+w_bus))/dx,1);

n_grid=n_grid_bus.*(1-core_reg_logical_ring)+n_grid_ring.*core_reg_logical_ring;
n_grid_bend=n_grid_bend_bus.*(1-core_reg_logical_ring)+n_grid_bend_ring.*core_reg_logical_ring;

% n_grid=max(n_grid_ring,n_grid_bus);
% n_grid_bend=max(n_grid_bend_ring,n_grid_bend_bus);

centering_shift=0.5*(w_bus+gap)/dx;
n_grid=circshift(n_grid,centering_shift,1);
n_grid_bend=circshift(n_grid_bend,centering_shift,1);

eps=n_grid.^2;
eps_bend=n_grid_bend.^2;

% Now we stretch out the mesh at the boundaries:
[x,y,xc,yc,dx_str,dy_str] = stretchmesh(x,y,[100,100,100,100],[4,4,4,4]);

%% Mode Calculation
guess=max(max(n_grid_bend));

boundary='0000'; % See wgmodes.m for details

t_start=tic();
[Hx,Hy,neff]=wgmodes(lambda0,guess,N_calc,dx_str,dy_str,eps_bend,boundary);

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
%     modes_Efields{i,1}=Efields_i/e_normfac;
%     modes_Efields{i,2}=Hfields_i/h_normfac;

    % Normalization method: \integ{E \times conj(H)}dA = 1
    normfac=sqrt(sum(sum(dotprod_2dbeams(mode_fields,conj_mode_struct(mode_fields),1)))*dx*dx);
    
    mode_fields.E=Efields_i/normfac;
    mode_fields.H=Hfields_i/normfac;
    
    modes_EHfields{i}=mode_fields;
end

clear ex ey ez hx hy hz Efields_i Ex Ey Ez Hfields_i Hx Hy Hz

[X,Y]=meshgrid(xc,flip(yc));
[Phi_grid,R_grid]=cart2pol(X,Y);

if plot_true
    close(gcf);
    chosenMode_Efields=modes_EHfields{plot_mode_ind}; chosenMode_Efields=chosenMode_Efields.E;
    plot_mode_Efields(X,Y,chosenMode_Efields,'jet',h_params,0);
    sgtitle(['Fundamental Hybrid Mode, Gap = ',num2str(gap),' nm']);
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
spat_arrays.n_grid_bend=n_grid_bend;
spat_arrays.h_params=h_params;
spat_arrays.hsliceBus=hsliceBus; spat_arrays.wsliceBus=wsliceBus;
spat_arrays.hsliceRing=hsliceRing; spat_arrays.wsliceRing=wsliceRing;