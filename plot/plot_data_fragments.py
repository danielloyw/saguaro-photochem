# -*- coding: utf-8 -*-
"""
Created on Sun Aug  9 13:01:14 2015

@author: rogeryelle
"""

sdir = ' '
read,sdir,prompt=' Enter run ID ... >'

; read list of nspecies

nmolec=read_molec_tchem1D('../wrk/nmolecules.dat')

imolec=read_molec_tchem1D('../wrk/imolecules.dat')

molec = [nmolec,imolec]

atm=read_atm_tchem1D('../runs/'+strtrim(sdir,2)+'/output/atm1D.out')

alt = 1.E-5*atm.alt
den = atm.den 
name = atm.name
nsp = n_elements(name)
nlev = n_elements(alt)
mole = fltarr(nlev,nsp)

for nm = 0, nsp-1 do begin
   mole[0:nlev-1,nm] = den[0:nlev-1,nm]/den[0:nlev-1,0]
endfor

print,name[0:nsp-1],format='(6(2X,A12,1X))'

name = strtrim(name,2)
neutr = n_elements(nmolec)
nsp = neutr

loadct,13
   set_plot, 'ps'
   mcharsz = 0.7
   ncharsz = 0.5
   textsize=0.5
   
   !p.thick=3.0
   !x.thick=3.0
   !y.thick=3.0
   !p.charthick=3.0

!p.multi=0
device,file='../runs/'+strtrim(sdir,2)+'/plots/model_data.ps',/color,/portrait
position=[0.15,0.15,0.90,0.90]


dum='    '
get_lun, lu
openr,lu,'../data/observations/mol_obs.txt'
nmol=FILE_LINES('../data/observations/mol_obs.txt')
mol_name=strarr(nmol)
obs=intarr(nmol)
FOR n=0,nmol-1 DO BEGIN
   readf,lu,dum
   mol_name(n)=strtrim(dum,2)
ENDFOR
close, lu
free_lun, lu

FOR k=0,nmol-1 DO BEGIN
   x=where(mol_name(k) eq name)
   IF (x ne -1) THEN BEGIN
      obs(k)=x
   ENDIF
   IF (x eq -1) THEN BEGIN
      print, mol_name(k), ' not found'
      obs(k)=-1
   ENDIF
ENDFOR

FOR l=0,nmol-1 DO BEGIN
   umin=1.e-1
   umax=1.e-12
   cmin=1.e-1
   cmax=1.e-12
   imin=1.e-1
   imax=1.e-12
   mmin=1.e-1
   mmax=1.e-12
   sfilename='../data/observations/'+strtrim(mol_name(l),2)+'_obs.txt'
   sline= ' '
   nc=0
   nu=0
   ni=0
   get_lun,lu
   openr,lu,sfilename
   readf,lu,sline, nc, FORMAT='(A5,I2)'
   readf,lu,sline, nu, FORMAT='(A5,I2)'
   readf,lu,sline, ni, FORMAT='(A5,I2)'
   IF (nc gt 0) THEN BEGIN
      lat_c=strarr(nc)
      ptrs_c=ptrarr(nc)
      type_c=intarr(nc)
      FOR n=0,nc-1 DO BEGIN
         nl=0
         readf,lu,sline
         readf,lu,sline
         mol=sline
         readf,lu,sline
         lat_c(n)=sline
         readf,lu,sline
         readf,lu,sline
         readf,lu,sline
         readf,lu,sline
         readf,lu,sline
         readf,lu,nl,dum1
         type_c(n)=dum1
         dummy=fltarr(6,nl)
         readf,lu,sline
         FOR m=0,nl-1 DO BEGIN
            readf, lu, dum1, dum2, dum3, dum4, dum5, dum6
            dummy(0,m)=dum1
            dummy(1,m)=dum2
            dummy(2,m)=dum3
            dummy(3,m)=dum4
            dummy(4,m)=dum5
            dummy(5,m)=dum6
         ENDFOR
         IF (min(dummy(3,*)) lt cmin) THEN BEGIN
            cmin=min(dummy(3,*))
         ENDIF
         IF (max(dummy(3,*)) gt cmax) THEN BEGIN
            cmax=max(dummy(3,*))
         ENDIF
         ptrs_c(n)=ptr_new(dummy)
      ENDFOR 
   ENDIF
   
   IF (nu gt 0) THEN BEGIN
      readf,lu,sline
      readf,lu,sline
      readf,lu,sline
      readf,lu,nl, dum
      uvis=fltarr(3,nl)
      FOR n=0,nl-1 DO BEGIN
         readf,lu, dum1, dum2, dum3
         uvis(0,n)=dum1
         uvis(1,n)=dum2
         uvis(2,n)=dum3
      ENDFOR
      umin=min(uvis(1,*))
      umax=max(uvis(1,*))
   ENDIF
   
   IF (ni gt 0) THEN BEGIN
      type_i=intarr(ni)
      ptrs_i=ptrarr(ni)
      FOR n=0,ni-1 DO BEGIN
         readf,lu,sline
         readf,lu,sline
         readf,lu,sline
         readf,lu,nl,dum
         type_i(n)=dum
         dummy1=fltarr(3,nl)
         FOR m=0,nl-1 DO BEGIN
            readf, lu, dum1, dum2, dum3
            dummy1(0,m)=dum1
            dummy1(1,m)=dum2
            dummy1(2,m)=dum3
         ENDFOR
         IF (min(dummy1(1,*)) lt imin) THEN BEGIN
            imin=min(dummy1(1,*))
         ENDIF
         IF (max(dummy1(1,*)) gt imax) THEN BEGIN
            imax=max(dummy1(1,*))
         ENDIF
         ptrs_i(n)=ptr_new(dummy1)
      ENDFOR
   ENDIF
   close, lu
   free_lun, lu

  IF (nc+ni+nu gt 0) THEN BEGIN  
   x=obs(l)
   If (x ge 0) THEN BEGIN
      mmin=min(mole(*,x))
      mmax=max(mole(*,x))
   ENDIF
   xmin_arr=[umin,imin,cmin,mmin]
   xmax_arr=[umax,imax,cmax,mmax]
   xmin=min(xmin_arr)/10
   xmax=max(xmax_arr)*10
   ;device, file=strtrim(mol_name(l),2)+'_mole.ps', /color
;plot, (*ptrs(0))(4,*), (*ptrs(0))(1,*), /xlog,/nodata, /ylog, yrange=[1.5,1e-3], xrange=[1e-9,5e-6], xstyle=1, $
;xtitle='Mixing Ratio', ytitle='Pressure (mbar)'
   IF (nc gt 0) THEN BEGIN
      plot, (*ptrs_c(0))(4,*), (*ptrs_c(0))(0,*), /xlog,/nodata, ystyle=1,yrange=[0,1200], xrange=[xmin,xmax], xstyle=1, $
            xtitle='Mixing Ratio', ytitle='Altitude (km)',title=strtrim(mol_name(l),2)
   ENDIF
   IF (nc eq 0) THEN BEGIN
      plot, (*ptrs_i(0))(1,*), (*ptrs_i(0))(0,*), /xlog,/nodata, ystyle=1,yrange=[0,1200], xrange=[xmin,xmax], xstyle=1, $
            xtitle='Mixing Ratio', ytitle='Altitude (km)',title=strtrim(mol_name(l),2)
   ENDIF
   plots, [xmin, xmax], [500,500], linestyle=2
   plots, [xmin, xmax], [1000,1000], linestyle=2
   xyouts, xmin*1.1, 460, 'CIRS'
   xyouts, xmin*1.1, 960, 'UVIS (-6S)'
   xyouts, xmin*1.1, 1010, 'INMS'
   IF (nc gt 0) THEN BEGIN
      FOR n=0,nc-1 DO BEGIN
         oplot, (*ptrs_c(n))(4,*), (*ptrs_c(n))(0,*), color=(300/nc)*(n)
         xyouts,0.15,n*0.03+0.15, lat_c(n), color=(300/nc)*(n),/normal
      ENDFOR
   ENDIF
   IF (nu gt 0) THEN BEGIN
      oplot, uvis(1,*), uvis(0,*), psym=-1
   ENDIF
   IF (ni gt 0) THEN BEGIN
      FOR n=0, ni-1 DO BEGIN
         IF (type_i(n) gt 0) THEN BEGIN
            oplot, (*ptrs_i(n))(1,*), (*ptrs_i(n))(0,*), psym=2,color=(300/ni)*(n)
         ENDIF
         IF (type_i(n) eq 0) THEN BEGIN   
            plotsym, 6,2, thick=3
            oplot, (*ptrs_i(n))(1,*), (*ptrs_i(n))(0,*), psym=8,color=(300/ni)*(n)
         ENDIF      
      ENDFOR
   ENDIF
; xyouts, 0.15, 0.9, strtrim(mol_name(l),2), /normal
x=obs(l)
If (x ge 0) THEN BEGIN
   oplot, mole(*,x), alt
ENDIF
   
ENDIF
ENDFOR   
device, /close_file  


end
