close all
clearvars
clc

% mydir=pwd;
% idcs=strfind(mydir,'/');
% % parentdir=mydir(1:idcs(end)-1);
%% Physical constants
addpath 'D:\Personal GDrive\Research stuff\Software and Packages\MATLAB Functions';
fundamental_constants;
hbar=physconsts.hbar; % in Joule seconds
epsilon0=physconsts.eps0;
c=physconsts.c;
mu0=1/(c^2*epsilon0);

%% Waveguide Parameters
lambda0=1550;

n1 = sellmeier('sio2',lambda0*10^-9);          % SiO2 lower cladding
n2 = sellmeier('sin',lambda0*10^-9);          % SiN core
n3 = sellmeier('hto',lambda0*10^-9);          % HTO upper cladding

% NOTE: All length dimensions here are in nm
h1 = 8000;           % SiO2 cladding height (referring to this as "height" assuming horizontal orientation of longer dimension of core)
h2 = 640;           % SiN core height
h3 = 8000;           % HTO cladding height
w = 1800;            % waveguide full-width
side = 2000;         % space on side of waveguide

dx=5;
dy=dx;
%% Index Grid

fprintf (1,'generating index mesh...\n');
[x,y,xc,yc,nx,ny,eps,edges] = ...
    waveguidemesh([n1,n2,n3],[h1,h2,h3],h2,w/2,side,dx,dy);
[x,y,xc,yc,dx_str,dy_str] = stretchmesh(x,y,[200,200,200,0],[1.5,1.5,1.5,1]);

% dx_str=dx; dy_str=dx;
n_grid=sqrt(eps);

%% Mode Calculation #1 - 000A (Ex symmetric about y-axis "A modes")

N_calc_AModes=2;
plot_mode_ind=1;
guess=max(max(n_grid));

boundary='000A'; % See wgmodes.m for details

t_start=tic();
[Hx,Hy,neff]=wgmodes(lambda0,guess,N_calc_AModes,dx_str,dy_str,eps,boundary);

neff_AModes=unique(abs(neff),'stable'); % Since we're not including scattering loss, this is okay
clear neff
keep_AModes_idx=find(neff_AModes>min([n1,n2,n3])); % Only keeping modes that would be guided under conventional TIR conditions
neff_AModes=neff_AModes(keep_AModes_idx);
N_AModes=length(neff_AModes);
toc(t_start);

AModes_fields=cell(N_AModes,2); % Column 1 corresponds to E fields, column 2 to H fields
for i=1:N_AModes
    ii=keep_AModes_idx(i);
    hx=Hx(:,:,ii);
    hy=Hy(:,:,ii);
    [hz,ex,ey,ez]=postprocess(lambda0,neff_AModes(ii),hx,hy,dx_str,dy_str,eps,boundary);
    
    % Cropping H fields to match the size of E fields
    hx=hx(1:end-1,1:end-1);
    hy=hy(1:end-1,1:end-1);
    hz=hz(1:end-1,1:end-1);
        
    ex_4q=[flipud(ex);ex];
    ey_4q=[-1*flipud(ey);ey];
    ez_4q=[flipud(ez);ez]; % Assuming here that Ez is symmetric about the y-axis...which it may or may not be
    
    hx_4q=[-1*flipud(hx);hx];
    hy_4q=[flipud(hy);hy];
    hz_4q=[flipud(hz);hz]; % Assuming here that Hz is symmetric about the y-axis...which it may or may not be
    
    Efields_i(:,:,1)=ex_4q;
    Efields_i(:,:,2)=ey_4q;
    Efields_i(:,:,3)=ez_4q;
    
    Hfields_i(:,:,1)=hx_4q;
    Hfields_i(:,:,2)=hy_4q;
    Hfields_i(:,:,3)=hz_4q;
    
    AModes_fields{i,1}=Efields_i;
    AModes_fields{i,2}=Hfields_i;
end

clear ex ey ez hx hy hz Efields_i Hfields_i Ex Ey Ez Hx Hy Hz

xc_4q=[flipud(-1*xc);xc];
yc_4q=yc;
[X,Y]=meshgrid(xc_4q,flip(yc_4q));
[Phi_grid,R_grid]=cart2pol(X,Y);

plot_mode_sixfields(X,Y,AModes_fields{plot_mode_ind,1},AModes_fields{plot_mode_ind,2},'turbo',h1);


%% Mode Calculation #2 - 000S (Ex anti-symmetric about y-axis "S modes")

N_calc_SModes=N_calc_AModes;
guess=max(max(n_grid));

boundary='000S'; % See wgmodes.m for details

t_start=tic();
[Hx,Hy,neff]=wgmodes(lambda0,guess,N_calc_SModes,dx_str,dy_str,eps,boundary);

neff_SModes=unique(abs(neff),'stable'); % Since we're not including scattering loss, this is okay
clear neff
keep_SModes_idx=find(neff_SModes>min([n1,n2,n3]));
neff_SModes=neff_SModes(keep_SModes_idx);
N_SModes=length(neff_SModes);
toc(t_start);

SModes_fields=cell(N_SModes,2); % Column 1 corresponds to E fields, column 2 to H fields
for i=1:N_SModes
    ii=keep_SModes_idx(i);
    hx=Hx(:,:,ii);
    hy=Hy(:,:,ii);
    [hz,ex,ey,ez]=postprocess(lambda0,neff_SModes(ii),hx,hy,dx_str,dy_str,eps,boundary);

    % Cropping H fields to match the size of E fields
    hx=hx(1:end-1,1:end-1);
    hy=hy(1:end-1,1:end-1);
    hz=hz(1:end-1,1:end-1);

    ex_4q=[-1*flipud(ex);ex];
    ey_4q=[flipud(ey);ey];
    ez_4q=[flipud(ez);ez]; % Assuming here that Ez is symmetric about the y-axis...which it may or may not be
    
    hx_4q=[flipud(hx);hx];
    hy_4q=[-1*flipud(hy);hy];
    hz_4q=[flipud(hz);hz]; % Assuming here that Hz is symmetric about the y-axis...which it may or may not be
    
    Efields_i(:,:,1)=ex_4q;
    Efields_i(:,:,2)=ey_4q;
    Efields_i(:,:,3)=ez_4q;
    
    Hfields_i(:,:,1)=hx_4q;
    Hfields_i(:,:,2)=hy_4q;
    Hfields_i(:,:,3)=hz_4q;
    
    SModes_fields{i,1}=Efields_i;
    SModes_fields{i,2}=Hfields_i;
end

clear ex ey ez hx hy hz Efields_i Hfields_i Ex Ey Ez Hx Hy Hz

xc_4q=[flipud(-1*xc);xc];
yc_4q=yc;
[X,Y]=meshgrid(xc_4q,flip(yc_4q));
[Phi_grid,R_grid]=cart2pol(X,Y);

plot_mode_sixfields(X,Y,SModes_fields{plot_mode_ind,1},SModes_fields{plot_mode_ind,2},'turbo',h1);

%% Combining field

N_tot=N_AModes+N_SModes;

vector_mode_fields_unsorted=[AModes_fields;SModes_fields];
vector_mode_fields=vector_mode_fields_unsorted;

mode_neff=[neff_AModes;neff_SModes];
[~,sortidx]=sort(mode_neff);

for i=1:length(mode_neff)
    vector_mode_fields{i}=vector_mode_fields_unsorted{sortidx(i),:};
end

fprintf(['Total number of modes supported in the waveguide = ',num2str(N_tot),'\n'])
%% Exporting
spat_arrays.X=X;
spat_arrays.Y=Y;
spat_arrays.R_grid=R_grid;
spat_arrays.Phi_grid=Phi_grid;
spat_arrays.xc_4q=xc_4q;
spat_arrays.dx=dx;
spat_arrays.eps=eps;
spat_arrays.n_grid=n_grid;

% save([parentdir,'/channelwg_1550nmlam_',num2str(w),'nm_by_',num2str(h2),'nm_unnorm'],'spat_arrays','vector_mode_fields','mode_neff');