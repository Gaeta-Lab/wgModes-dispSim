function [neff,varargout]=modecalc_fullgrid_step_func2(h_params,lambda0,N_calc,w)
%% Index Grid

if nargin==3
    modefield_true=0;
end

lowCladMat=h_params.lowCladMat;
coreMat=h_params.coreMat;
upperCladMat=h_params.upperCladMat;

h1=h_params.h1; h2=h_params.h2; h3=h_params.h3; side=h_params.side; %w=h_params.w;
dx=h_params.dx; dy=h_params.dy;
bend_radius_nm=h_params.bend_radius_nm;

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
guess=max(max(n_grid_bend));

boundary='0000'; % See wgmodes.m for details

% t_start=tic();
[Hx,Hy,neff]=wgmodes(lambda0,guess,N_calc,dx_str,dy_str,eps_bend,boundary);

neff=unique(abs(neff),'stable'); % Since we're not including scattering loss, this is okay
keep_modes_idx=find(neff>min([n1,n2,n3])); % Only keeping modes that would be guided under conventional TIR conditions
neff=neff(keep_modes_idx);
N_modes=length(neff);

if N_modes<N_calc
    neff=[neff;zeros(N_calc-N_modes,1)];
end

modes_EHfields=cell(N_modes,1);
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

        % Normalization method: \integ{E \times conj(H)}dA = 1
        normfac=sqrt(sum(sum(dotprod_2dbeams(mode_fields,conj_mode_struct(mode_fields),1)))*dx*dx);

        mode_fields.E=Efields_i/normfac;
        mode_fields.H=Hfields_i/normfac;

        modes_EHfields{i}=mode_fields;
    end
    
[X,Y]=meshgrid(xc,flip(yc));

clear ex ey ez hx hy hz Efields_i Ex Ey Ez Hfields_i Hx Hy Hz


h_params.w=w;

spat_arrays.X=X; spat_arrays.Y=Y;
spat_arrays.xc=xc; spat_arrays.dx=dx;
spat_arrays.eps=eps; spat_arrays.n_grid=n_grid; spat_arrays.n_grid_bend=n_grid_bend;
spat_arrays.core_reg_logical=core_reg_logical; spat_arrays.h_params=h_params;
spat_arrays.N_modes=N_modes;

varargout{1}=modes_EHfields;
varargout{2}=spat_arrays;
