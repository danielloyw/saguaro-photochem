  fid = fopen('eimpact.dat','r');
  
  neng = fscanf(fid,'%d',1);
  nabs_thk = fscanf(fid,'%d',1);
  nabs_thn = fscanf(fid,'%d',1);
  nbrnmax = fscanf(fid,'%d',1);
  nprdmax = fscanf(fid,'%d\n',1);
  nabs = nabs_thk + nabs_thn;
  fgetl(fid);
  eng = fscanf(fid,'%e\n',neng);
  fgetl(fid); 
  deng = fscanf(fid,'%e\n',neng);
 
  for na = 1:1  
     xname=fscanf(fid,'%s',1);
     nexcite=fscanf(fid,'%d',1);
     ndissoc=fscanf(fid,'%d',1);
     nioniz=fscanf(fid,'%d\n',1);
     nbrnch = ndissoc + nioniz;
     ncs = nexcite + ndissoc + nioniz;
     fgetl(fid);     
     crs_tot = fscanf(fid,'%e\n',neng);
     fgetl(fid);     
     crs_inel = fscanf(fid,'%e\n',neng);
     ecs = zeros(neng,ncs);     
     nx = 0;
     for ni = 1:nexcite
        nx = nx + 1;
        yname = fscanf(fid,'%s',1);
        yenrg = fscanf(fid,'%f',1);
        sdum = fscanf(fid,'%s\n',1);
        ecs(1:neng,nx) = fscanf(fid,'%e\n',neng);
     end
     for ni = 1:nbrnch
        nx = nx + 1;
        yname = fscanf(fid,'%s',1);
        yenrg = fscanf(fid,'%f',1);
        sdum = fscanf(fid,'%s\n',1);
        fgetl(fid);
        ecs(1:neng,nx) = fscanf(fid,'%e\n',neng);
     end
  end
  
fclose(fid);

figure; 
hold on; 
grid on; 
set(gca,'fontsize',14,'xlim',[0.1,1000],'ylim',[1.E-20,1.E-13],'XMinorTick','on',...
    'YMinorTick','on','FontWeight','bold','xscale','log','yscale','log'); 
xlabel('Energy (eV)','FontSize',14,'FontWeight','bold'); 
ylabel('Cross Section (cm^2)','FontSize',14,'FontWeight','bold');

plot(gca,eng,crs_tot,'-r','linewidth',2);
plot(gca,eng,crs_inel,'-b','linewidth',2);
for nx=1:ncs
plot(gca,eng,ecs(1:neng,nx),'-k','linewidth',2);
end

saveas(gcf,'elctrn.fig');
