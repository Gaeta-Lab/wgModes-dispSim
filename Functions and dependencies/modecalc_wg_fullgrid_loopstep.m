%% Index Grid

n1 = sellmeier(lowCladMat,lambda0*10^-9);          % SiO2 lower cladding
n2 = sellmeier(coreMat,lambda0*10^-9);          % SiN core
n3 = sellmeier(upperCladMat,lambda0*10^-9);          % HTO upper cladding

% fprintf (1,'generating index mesh...\n');
[x,y,xc,yc,nx,ny,eps,edges] = ...
    waveguidemeshfull([n1,n2,n3],[h1,h2,h3],h2,w/2,side,dx,dy); 

% Now we stretch out the mesh at the boundaries:
[x,y,xc,yc,dx_str,dy_str] = stretchmesh(x,y,[100,100,100,100],[4,4,4,4]);
n_grid=sqrt(eps);

% "Simulating" bending by adding a slope to the core index profile
core_reg_logical=(n_grid>0.9*n2);
bend_pad_width=400; % Width of region in cladding on either side of core that I simulate bending for
bend_pad_ind=round(bend_pad_width/dx);

[co,ro,~]=find(core_reg_logical==1);
hslice=unique(ro);
wslice=unique(co');
wslice=wslice(1)-bend_pad_ind:wslice(end)+bend_pad_ind;

n_grid_bend=n_grid;

bend_factor=(1+(xc(wslice)-xc(wslice(round(length(wslice)/2))))/bend_radius_nm);
bend_factor=bend_factor';

bend_matrix=1+0*hslice*wslice;
bend_matrix=bend_matrix.*bend_factor; bend_matrix=bend_matrix';

n_grid_bend(wslice,hslice)=n_grid(wslice,hslice).*bend_matrix;
eps_bend=n_grid_bend.^2;

%% Mode Calculation

plot_mode_ind=1;
guess=max(max(n_grid_bend));

boundary='0000'; % See wgmodes.m for details

% t_start=tic();
[Hx,Hy,neff]=wgmodes(lambda0,guess,N_calc,dx_str,dy_str,eps_bend,boundary);

neff=unique(abs(neff),'stable'); % Since we're not including scattering loss, this is okay
keep_modes_idx=find(neff>min([n1,n2,n3])); % Only keeping modes that would be guided under conventional TIR conditions
neff=neff(keep_modes_idx);
N_modes=length(neff);
% toc(t_start);
