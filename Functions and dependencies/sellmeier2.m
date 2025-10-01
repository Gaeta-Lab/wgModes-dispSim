function out = sellmeier2(mat,wlength)
eV_to_nm = 1239.842; % nm*eV
wl_nm = wlength*1e9;
wl_um = wlength*1e6;

% Si: LPCVD - slightly wrong
if strcmpi(mat,'sin')
    out = sqrt(1+253.86./(9.1609.^2-(eV_to_nm./wl_nm).^2)+0.040314./(0.001.^2-(eV_to_nm./wl_nm).^2));

    %sqrt(1+253.86/(9.1609^2-(1239.842/wl_nm)^2)+0.040314/(0.001^2-(1239.842/wl_nm)^2));

% SiO2: Thermal SiO2
elseif strcmpi(mat,'sio2')
    out = sqrt(1.4124+85.493./(11.061.^2-(eV_to_nm./wl_nm).^2)+0.015888./(0.10237.^2-(eV_to_nm./wl_nm).^2));

% Si:
elseif strcmpi(mat,'si')
    out = sqrt(1+10.6684293/(1-(0.301516485/wl_um).^2)+0.0030434748/(1-(1.13475115/wl_um).^2)+1.54133408/(1-(1104/wl_um).^2));

% HTO SiO2: (its own process)
elseif strcmpi(mat,'hto')
    % out = sqrt(1+0.31197/(1-(5836.94547 [nm]/WL)^2)+0.056473676/(1-(735223.771[nm]/WL)^2)+1.16221658/(1-(100.475642[nm]/WL)^2));
    out = sqrt(1+0.31197./(1-(5836.94547./wl_nm).^2)+0.056473676./(1-(735223.771./wl_nm).^2)+1.16221658./(1-(100.475642./wl_nm).^2));
    
% MgF2
elseif strcmpi(mat,'mgf2')
    out = sqrt(1+0.48755108/(1-(0.04338408/wl_um)^2)+0.39875031/(1-(0.09461442/wl_um)^2)+2.3120353/(1-(23.793604/wl_um)^2));
    % Diamond
elseif strcmpi(mat,'diamond')
    out = sqrt(1+0.3306/(1-(0.1750/wl_um)^2)+4.3356/(1-(0.1060/wl_um)^2));

    % SiN: PECVD after anneal
elseif strcmpi(mat,'sin_pe_be')
    % Data from following paper:
    % https://onlinelibrary.wiley.com/doi/full/10.1002/lpor.202200544
    out = sqrt( 1 + 2.61.*wl_nm.^2./(wl_nm.^2 -139.77.^2) +...
        1.11e9.*wl_nm.^2./(wl_nm.^2 -(2.51e8).^2) );

    % SiN: PECVD after anneal
elseif strcmpi(mat,'sin_pe')
    % Data from following paper:
    % https://onlinelibrary.wiley.com/doi/full/10.1002/lpor.202200544
    out = sqrt( 1 + 2.97.*wl_nm.^2./(wl_nm.^2 -144.86.^2) +...
        1.57e9.*wl_nm.^2./(wl_nm.^2 -(3.8e8).^2) );

elseif strcmpi(mat,'sio2_pe')
    % SiO2: PECVD
    out = sqrt(1+0.72158.*wl_um.^2/(wl_um.^2-(0.089746).^2)+...
        0.43225.*wl_um.^2./(wl_um.^2-(0.099851).^2)+...
        0.66606.*wl_um.^2/(wl_um.^2-(9.7123)^2));
elseif strcmpi(mat,'aln')
    % AIM AlN
    out = sqrt( 1 + 2.1070./(1 - (99.1702/wl_nm).^2 ) + 1.0428./(1 - (170.1432/wl_nm).^2 )  );
elseif strcmpi(mat,'alo')
    a1=1.761e+10; a2=1.561e+04; a3=1.724e+10;
    b1=-2.093e+10; b2=2.116e+04; b3=-2.057e+10;
    out=sqrt(1+a1/(wl_nm^2-b1)+a2/(wl_nm^2-b2)+a3/(wl_nm^2-b3));
    % Diamond
% Correct 2015 Luke SiN
elseif strcmpi(mat,'ksin')
    A = 3.0249;
    B = 40314;
    %B = 40341;
    X = 0.1353406;
    Y = 1347.3;
    out = sqrt( 1 + (A.*wl_um.^2)./(wl_um.^2 - X.^2) + ...
        (B.*wl_um.^2)./(wl_um.^2 - Y.^2) );
% Raw 2015 Luke SiN Spline
elseif strcmpi(mat,'ssin')
    wl_um2 = wl_um;
    load('nsin.mat')
    out = spline(wl_um,n,wl_um2);
% 2025 SiN
elseif strcmpi(mat,'nsin')
    A = 2.948;
    X = 0.14;
    B = 1.7669e7;
    Y = 1e8;
    out = sqrt(1 + (A*wl_um.^2)./(wl_um.^2-X.^2) + ...
        (B*wl_um.^2)./(wl_um.^2 - Y.^2)   );
% Raw 2025 SiN Spline
elseif strcmpi(mat,'ssin2')

    data = readmatrix('B4_sin365_anneal_fit_wLoss.txt');

    lam = data(:,1);
    n = data(:,2);
    out = spline(lam,n,wl_nm);
% 2025 HTO
elseif strcmpi(mat,'hto2')
    A = 0.733;
    Einf = 1.372;
    B = 0.11107;
    E = 0.01038;
    out = sqrt(Einf + A.*wl_um.^2./(wl_um.^2-B.^2) - E.*wl_um.^2);

end

%{
function out = sellmeier(mat,wlength)
eV_to_nm = 1239.842; % nm*eV
wl_nm = wlength*1e9;
wl_um = wlength*1e6;

if strcmpi(mat,'sin')
    % For SiN:
    out = sqrt(1+253.86/(9.1609^2-(eV_to_nm/wl_nm)^2)+0.040314/(0.001^2-(eV_to_nm/wl_nm)^2));
elseif strcmpi(mat,'sio2')
    % For SiO2:
    out = sqrt(1.4124+85.493/(11.061^2-(eV_to_nm/wl_nm)^2)+0.015888/(0.10237^2-(eV_to_nm/wl_nm)^2));
elseif strcmpi(mat,'si')
    % For Si:
    out = sqrt(1+10.6684293/(1-(0.301516485/wl_um).^2)+0.0030434748/(1-(1.13475115/wl_um).^2)+1.54133408/(1-(1104/wl_um).^2));
elseif strcmpi(mat,'hto')
    % For HTO SiO2:
    % out = sqrt(1+0.31197/(1-(5836.94547 [nm]/WL)^2)+0.056473676/(1-(735223.771[nm]/WL)^2)+1.16221658/(1-(100.475642[nm]/WL)^2));
    out = sqrt(1+0.31197/(1-(5836.94547/wl_nm).^2)+0.056473676/(1-(735223.771/wl_nm).^2)+1.16221658/(1-(100.475642/wl_nm)^2));
    % MgF2
elseif strcmpi(mat,'mgf2')
    out = sqrt(1+0.48755108/(1-(0.04338408/wl_um)^2)+0.39875031/(1-(0.09461442/wl_um)^2)+2.3120353/(1-(23.793604/wl_um)^2));
    % AlN
elseif strcmpi(mat,'aln')
    a1=8.164e+08; a2=8.205e+08; a3=4.656e+04;
    b1=-5.2e+08; b2=-5.15e+08; b3=2.653e+04;
    out=sqrt(1+a1/(wl_nm^2-b1)+a2/(wl_nm^2-b2)+a3/(wl_nm^2-b3));
    % AlO
elseif strcmpi(mat,'alo')
    a1=1.761e+10; a2=1.561e+04; a3=1.724e+10;
    b1=-2.093e+10; b2=2.116e+04; b3=-2.057e+10;
    out=sqrt(1+a1/(wl_nm^2-b1)+a2/(wl_nm^2-b2)+a3/(wl_nm^2-b3));
    % Diamond
elseif strcmpi(mat,'diamond')
    out = sqrt(1+0.3306/(1-(0.1750/wl_um)^2)+4.3356/(1-(0.1060/wl_um)^2));
end
%}