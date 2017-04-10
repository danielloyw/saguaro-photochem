q0 = 6.512E-14;

%
%  Read energy grid
%

egrid = csvread('egrid.csv');
neng = size(egrid,1);
eng = egrid(1:neng,1);

%
%  Read data from Itikawa
%

pv1 = csvread('CO2_V1.csv',1,0);
ndat = size(pv1,1);
edat = pv1(1:ndat,1);
cdat = 1.E-16*pv1(1:ndat,2);

%
%  Read fit parameters from Bhardwaj & Jain
%

params_tbl4 = csvread('Bhardwaj_table4.csv',1,1,[1,1,19,7]);
W = params_tbl4(1,1);
alpha = params_tbl4(1,2);
beta = params_tbl4(1,3);
WJ = params_tbl4(1,4);
Omega = params_tbl4(1,5);
F = params_tbl4(1,6);
AF = params_tbl4(1,7);


crs_v1 = zeros(neng,1);
indx = find((eng > edat(1)) & (eng <= edat(ndat)));
crs_v1(indx) = exp(interp1(edat,log(cdat),eng(indx)));
indx = find(eng>edat(ndat));
crs_v1(indx) = (q0*F/W^2)*((1.0 - (W./eng(indx)).^alpha).^beta).*(W./eng(indx)).^Omega;
