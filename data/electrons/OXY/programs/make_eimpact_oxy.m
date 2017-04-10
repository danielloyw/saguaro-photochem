q0 = 6.512E-14;

%
%  Read energy grid
%

egrid = csvread('egrid.csv');
neng = size(egrid,1);
eng = egrid(1:neng,1);
deng = egrid(1:neng,2);

%%
%  Total Elastic Cross Section
%

data = csvread('Total_Elastic.csv',1,0);
ndat = size(data,1);

crs_tot_elastic = zeros(neng,1);
crs_tot_elastic(1:5) = interp1(data(1:ndat,1),data(1:ndat,2),eng(1:5),'linear','extrap');
crs_tot_elastic(6:neng) = exp(interp1(data(1:ndat,1),log(data(1:ndat,2)),eng(6:neng)));


%%  Excitation Cross Sections

nexcit=19;
crs_excit = zeros(neng,19);
ethreshold = zeros(19,1);

%% Vibrational 

%
%  Read fit parameters from Bhardwaj & Jain
%

params_tbl4 = csvread('Bhardwaj_table4.csv',1,1,[1,1,19,7]);

%
%  Read V1 data from Itikawa
%

pv1 = csvread('CO2_V1.csv',1,0);
nd1 = size(pv1,1);
ed1 = pv1(1:nd1,1);
cd1 = 1.E-16*pv1(1:nd1,2);

W = params_tbl4(1,1);
alpha = params_tbl4(1,2);
beta = params_tbl4(1,3);
WJ = params_tbl4(1,4);
Omega = params_tbl4(1,5);
F = params_tbl4(1,6);

indx = find((eng > ed1(1)) & (eng <= ed1(nd1)));
crs_excit(indx,1) = exp(interp1(ed1,log(cd1),eng(indx)));
ind1 = find(eng>ed1(nd1));
crs_excit(ind1,1) = (q0*F/W^2)*((1.0 - (W./eng(ind1)).^alpha).^beta).*(W./eng(ind1)).^Omega;
ethreshold(1) = W;

%
%  Read V2 data from Itikawa
%

pv2 = csvread('CO2_V2.csv',1,0);
nd2 = size(pv2,1);
ed2 = pv2(1:nd2,1);
cd2 = 1.E-16*pv2(1:nd2,2);

W = params_tbl4(2,1);
alpha = params_tbl4(2,2);
beta = params_tbl4(2,3);
WJ = params_tbl4(2,4);
Omega = params_tbl4(2,5);
F = params_tbl4(2,6);

indx = find((eng > ed2(1)) & (eng <= ed2(nd2)));
crs_excit(indx,2) = exp(interp1(ed2,log(cd2),eng(indx)));
ind2 = find(eng>ed2(nd2));
crs_excit(ind2,2) = (q0*F/W^2)*((1.0 - (W./eng(ind2)).^alpha).^beta).*(W./eng(ind2)).^Omega;
ethreshold(2) = W;

%
%  Read V3 data from Itikawa
%

pv3 = csvread('CO2_V3.csv',1,0);
nd3 = size(pv3,1);
ed3 = pv3(1:nd3,1);
cd3 = 1.E-16*pv3(1:nd3,2);

W = params_tbl4(3,1);
alpha = params_tbl4(3,2);
beta = params_tbl4(3,3);
WJ = params_tbl4(3,4);
Omega = params_tbl4(3,5);
F = params_tbl4(3,6);

indx = find((eng > ed3(1)) & (eng <= ed3(nd3)));
crs_excit(indx,3) = exp(interp1(ed3,log(cd3),eng(indx)));
ind3 = find(eng>ed3(nd3));
crs_excit(ind3,3) = (q0*F/W^2)*((1.0 - (W./eng(ind3)).^alpha).^beta).*(W./eng(ind3)).^Omega;
ethreshold(3) = W;

%%
%  Excitation of electronic states
%

params_tbl4 = csvread('Bhardwaj_table4.csv',1,1,[1,1,19,7]);
for ns=4:19
    W = params_tbl4(ns,1);
    alpha = params_tbl4(ns,2);
    beta = params_tbl4(ns,3);
    WJ = params_tbl4(ns,4);
    Omega = params_tbl4(ns,5);
    F = params_tbl4(ns,6);
    AF = params_tbl4(ns,7);
    indx = find(eng>W);
    crs_excit(indx,ns) = (q0*F/W^2)*((1.0 - (W./eng(indx)).^alpha).^beta) ...
        .*(W./eng(indx)).^Omega;
    ethreshold(ns) = W;
end


%%
%  Ionization and Dissociative Ionization Cross Sections
%

crs_ioniz = zeros(neng,10);

%
%  e + CO2 --> CO2+(A), CO2+(B) from tables (fits in Bhardwaj & Jain are
%  poor)
%

pinz = csvread('Itikawa_table12.csv',1,0);
ndat = size(pinz,1);
edat = pinz(1:ndat,1);
indx = find(eng > edat(1));
crs_ioniz(indx,2) = exp(interp1(edat,log(1.E-16*pinz(1:ndat,2)),eng(indx)));
crs_ioniz(indx,3) = exp(interp1(edat,log(1.E-16*pinz(1:ndat,3)),eng(indx)));

%
%  Other channels from Bhardwaj & Jain
%

params_tbl3 = csvread('Bhardwaj_table3.csv',1,1,[1,1,14,11]);
[nstates_ioniz, nparams_tbl3] = size(params_tbl3);

ilist = [2,5,6,7,8,9,10,11];
crs_ioniz = zeros(neng,10);
ithreshold = zeros(10,1);
for nx=1:numel(ilist)
    ns = ilist(nx);
    W = params_tbl3(ns,1);
    I = params_tbl3(ns,2);
    K = params_tbl3(ns,3);
    KB = params_tbl3(ns,4);
    J = params_tbl3(ns,5);
    JB = params_tbl3(ns,6);
    JC = 0.0;
    GS = params_tbl3(ns,7);
    GB = params_tbl3(ns,8);
    TS = params_tbl3(ns,9);
    TA = params_tbl3(ns,10);
    TB = params_tbl3(ns,11);
    indx = find(eng>W);
    A = 1.0E-16*(K./(KB+eng(indx))).*log(eng(indx)./J + JB + JC./eng(indx));
    G = GS*(eng(indx)./(eng(indx)+GB));
    T0 = TS - (TA./(eng(indx)+TB));
    TM = 0.5*(eng(indx)-I);
    crs_ioniz(indx,ns-1) = A.*G.*(atan((TM-T0)./G)+atan(T0./G));    
    ithreshold(ns-1) = W;
end

crs_total_inelastic = zeros(neng,1);
for n=1:neng
    crs_total_inelastic(n)=crs_total_inelastic(n)+sum(crs_excit(n,1:19)) ...
        +sum(crs_ioniz(n,1:10));
end

sexcit=['VIB(010)    ';
        'VIB(100)    ';
        'VIB(001)    ';
        'ELCTRN_8.6  ';
        'ELCTRN_9.3  ';
        'ELCTRN_11.1 ';
        'ELCTRN_12.4 ';
        'ELCTRN_13.6 ';
        'ELCTRN_15.5 ';
        'ELCTRN_16.3 ';
        'ELCTRN_17.0 ';
        'ELCTRN_18.7 ';
        'OI_1304     ';
        'OI_1356     ';
        'CI_1279     ';
        'CI_1329     ';
        'CI_1561     ';
        'CI_1657     ';
        'COP(1N)     '];
    
sioniz=['CO2P        ';
        'CO2P        ';
        'CO2P        ';
        'CO2P        ';
        'COP         ';
        'OP          ';
        'CP          ';
        'CPP         ';
        'OPP         ';
        'CO2PP       '];

stitle=['CO2          + E            = CO2P         + E            +              +                ';...
        'CO2          + E            = CO2P         + E            +              +                ';...
        'CO2          + E            = CO2P         + E            +              +                ';...
        'CO2          + E            = CO2P         + E            +              +                ';...
        'CO2          + E            = COP          + O            + E            +                ';...
        'CO2          + E            = OP           + CO           + E            +                ';...
        'CO2          + E            = CP           + O2           + E            +                ';...
        'CO2          + E            = CPP          + O2           + E            +                ';...
        'CO2          + E            = OPP          + CO           + E            +                ';...
        'CO2          + E            = CO2PP        + E            +              +                '];
    
%%    WRITE EIMPACT.DAT
nexcit = 19;
ndissoc = 0;
nioniz = 10;
nbrnch = ndissoc + nioniz;
ncs = nexcit + ndissoc + nioniz;

fid = fopen('eimpact_O.dat','w');

    fprintf(fid,'%6d ',neng);
    nabs_thk = 1;
    fprintf(fid,'%6d',nabs_thk);
    nabs_thn = 0;
    fprintf(fid,'%6d',nabs_thn);
    nbrnmax = 0;
    fprintf(fid,'%6d',nbrnmax);
    nprdmax = 0;
    fprintf(fid,'%6d\n',nprdmax);
    fprintf(fid,'MEAN ENERGY IN BINS\n');
    fprintf(fid,' %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n',eng);
    fprintf(fid,'WIDTH OF ENERGY BINS\n'); 
    fprintf(fid,' %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n',deng);
    fprintf(fid,'CO2\n');
    fprintf(fid,'%6d %6d %6d\n',nexcit,ndissoc,nioniz);
    fprintf(fid,'TOTAL ELASTIC CROSS SECTION\n');
    fprintf(fid,' %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n',crs_tot_elastic);    
    fprintf(fid,'TOTAL INELASTIC CROSS SECTION\n');
    fprintf(fid,' %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n',crs_total_inelastic);
    for n = 1:nexcit
        fprintf(fid,'%12s %8.3f\n',sexcit(n,1:12),ethreshold(n));
        fprintf(fid,' %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n',crs_excit(1:neng,n));
    end
    for n = 1:nioniz
        fprintf(fid,'%12s %8.3f\n',sioniz(n,1:12),ithreshold(n));
        fprintf(fid,'%90s\n',stitle(n,1:90));
        fprintf(fid,' %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n',crs_ioniz(1:neng,n));
    end

fclose(fid);


