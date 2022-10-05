function dwd_scripts_4
%aus: /ocean-data/grashorn/nbs_fine_skag_juli_2014
%dwd2sflux.m
%Purpose: convert DWD wind (nc) to sflux for SELFE
%Authors: Johannes Pein (johannes.pein@uni-oldenburg.de); Joseph Zhang;
%slight modifications not affecting the output by Benjamin Jacob
%Date: Nov 2012
%Notes: does not work across years.
%Inputs: DWD nc files, constants specified near the beginning of this script.
%Outputs: wind/*.nc (sflux format)
% Last change 16.02.2021: make cloud shadow effectiveness a potential tuning parameter
% Last change 31.03.2015: replace varID calls by fixed number (which changed after gridcorrection)
                         % to varID call inquired by name                  
						 
clear all; close all;

month_day=[31 28 31 30 31 30 31 31 30 31 30 31];
year=2018; %2012 = Schaltjahr
start_mon=1; %starting month for model
start_day=2;

startdate=datenum(year,start_mon,start_day);
enddate=datenum(2018,12,31);

%ndays=88; %total # of run days
ndays=enddate-startdate+1; %total # of run days
startdate=startdate-1;
%ndays=1
% file settings
dwddir='/gpfs/work/ksddata/model/atmosphere/dwd/';
outdir='/gpfs/work/jacobb/data/RUNS/NorthSea/NorthSeaV2/sflux/';

% effective cloud shading effect
use_altered_cloud_eff=true 
cloud_shadow_eff_tuned=0.5
if use_altered_cloud_eff
	cloud_shadow_eff=cloud_shadow_eff_tuned;
else
	cloud_shadow_eff=0.62;
end
%dwddir='/home/benjamin/Desktop/Wind_Experiments/corrected_grid_winddata/dwd';
%outdir='/home/benjamin/Desktop/Wind_Experiments/corrected_grid_winddata/sflux';

%dwddir='/home/benjamin/Desktop/SCHISM/Forcing/windData/originalDWD'

if(rem(year,4)==0)
  month_day(2)=29;
  dschalt=1
else
    dschalt=2
end

mont=start_mon; day=start_day-1;
for d=1:ndays
  day=day+1;
  if(day>month_day(mont))
    day=day-month_day(mont);
    mont=mont+1;
    
    %reevaluate yeare
    yearstr=datestr(startdate+d,'YYYYmmdd');
    year=str2double(yearstr(1:4));
    %mont=str2double(yearstr(5:6)); %
    %day=str2double(yearstr(7:8));
    
  end
  
  
  %day=62
      

  mon_char=sprintf('%2.2d',mont);
  day_char=sprintf('%2.2d',day);
  if mont>2 %Schaltjahr solar radiation
    yday=datenum(year,mont,day,0,0,0)-datenum(year,1,1,0,0,0)+dschalt; %Julian day (Jan 1 is day 1)
  else
    yday=datenum(year,mont,day,0,0,0)-datenum(year,1,1,0,0,0)+1; %+1; %Julian day (Jan 1 is day 1)
  end
  datestr(yday)
  

  %DWD dir
  %fname=['/h/tmp/ksd/dwd' num2str(year) mon_char day_char '.nc']
  fname=fullfile(dwddir,['dwd' num2str(year) mon_char day_char '.nc']);
  
  %fname=['NC/dwd' num2str(year) mon_char day_char '.nc'];
  ncid1=netcdf.open(fname,'NC_NOWRITE');
  disp(['Reading ' fname]);
  
  
  time=netcdf.getVar(ncid1,netcdf.inqVarID(ncid1,'time')); %int32!!
  lon=netcdf.getVar(ncid1,netcdf.inqVarID(ncid1,'lon'));
  lat=netcdf.getVar(ncid1,netcdf.inqVarID(ncid1,'lat'));
  press=netcdf.getVar(ncid1,netcdf.inqVarID(ncid1,'slp'));
  u10=netcdf.getVar(ncid1,netcdf.inqVarID(ncid1,'u10'));
  v10=netcdf.getVar(ncid1,netcdf.inqVarID(ncid1,'v10'));
  t2=netcdf.getVar(ncid1,netcdf.inqVarID(ncid1,'t2')); %air temp. (K)
  dev2=netcdf.getVar(ncid1,netcdf.inqVarID(ncid1,'dev2')); %Dew point temp. (K)
  tcc=netcdf.getVar(ncid1,netcdf.inqVarID(ncid1,'tcc')); %Total cloud cover [0-1]
  prc=netcdf.getVar(ncid1,netcdf.inqVarID(ncid1,'precip')); %total precip. [kg/m^2]
  netcdf.close(ncid1);

  %Need to get new prc files
  %if(length(find(prc<0)) ~=0); error('prc<0'); end;
  ntime=length(time);
  dt=double(time(2)-time(1)); %in sec
  time2=0.:dt:(ntime-1.)*dt; %output time stamps
  
 keyboard

  %[X,Y] = meshgrid(lon,lat);
  %Be very careful about the orientation of quad; more later
  nx=length(lon); ny=length(lat);
  for i=1:nx
    xgrid(i,1:ny)=lon(i);
  end %for
  for j=1:ny
    ygrid(1:nx,j)=lat(j);
  end %for

  %Compute specific humidity (from Joanna Staneva)
  a1=6.107799961;
  a2=4.436518521e-1;
  a3=1.428945805e-2;
  a4=2.650648471e-4;
  a5=3.031240396e-6;
  a6=2.034080948e-8;
  a7=6.136820929e-11;
  const06=0.62198;
  ea= a1+a2*(dev2-273).*(a2+(dev2-273).*(a3+(dev2-273).*(a4+(dev2-273).* ...
      (a5+(dev2-273).*(a6+(dev2-273).*a7)))));
  ea=ea*100.0; % Conversion millibar --> Pascal
  spfh= const06*ea./(press-0.377*ea);
  indx=find(spfh<0);
  if(length(indx) ~=0); spfh(indx)=0; end; %error('Negative spfh'); end;

  %Solar rad
  for i=1:ntime
    %Had to conserve memory
    for j=1:ny
      hour=(i-1)*dt/3600;
      solar(:,j,i)=short_wave_radiation(yday,hour,xgrid(:,j),ygrid(:,j),tcc(:,j,i),cloud_shadow_eff);
    end %for j
  end %for i
  if(sum(sum(sum(isnan(solar)))) ~=0); error('nan in solar'); end;
  if(length(find(solar<0)) ~=0); error('Solar<0'); end;
 
  %_Downward_ longwave radiation (Josey et al. 2003a)
  bolz=5.6704e-8; %Stefan-Boltzmann const [W/m^2/K^4]
  longwv=bolz*(t2+10.77*tcc.^2+2.34*tcc-18.44).^4;

  %Now flip all matrices left-right as lat is from top-down; sflux code expect triangle area
  %(i,j), (i+1,j), and (i+1,j+1) is positive
  xgrid=fliplr(xgrid);
  ygrid=fliplr(ygrid);
  for i=1:ntime
    press(:,:,i)=fliplr(press(:,:,i));
    u10(:,:,i)=fliplr(u10(:,:,i));
    v10(:,:,i)=fliplr(v10(:,:,i));
    t2(:,:,i)=fliplr(t2(:,:,i));
    spfh(:,:,i)=fliplr(spfh(:,:,i));
    longwv(:,:,i)=fliplr(longwv(:,:,i));
    solar(:,:,i)=fliplr(solar(:,:,i));
    prc(:,:,i)=fliplr(prc(:,:,i));
  end %for i

  %Write *_air.nc and *_rad.nc
  %air
  ncid2=netcdf.create(fullfile(outdir,['air_' num2str(year) mon_char day_char '.nc']),'NC_CLOBBER');
  ncid3=netcdf.create(fullfile(outdir,['rad_' num2str(year) mon_char day_char '.nc']),'NC_CLOBBER');
  ncid4=netcdf.create(fullfile(outdir,['prc_' num2str(year) mon_char day_char '.nc']),'NC_CLOBBER');
  lon_dim=netcdf.defDim(ncid2,'nx_grid',nx);
  lat_dim=netcdf.defDim(ncid2,'ny_grid',ny);
  time_dim=netcdf.defDim(ncid2,'time',ntime);
  time_var=netcdf.defVar(ncid2,'time','float',time_dim);
  lon_var=netcdf.defVar(ncid2,'lon','float',[lon_dim lat_dim]);
  lat_var=netcdf.defVar(ncid2,'lat','float',[lon_dim lat_dim]);
  press_var=netcdf.defVar(ncid2,'prmsl','float',[lon_dim lat_dim time_dim]);
  u10_var=netcdf.defVar(ncid2,'uwind','float',[lon_dim lat_dim time_dim]);
  v10_var=netcdf.defVar(ncid2,'vwind','float',[lon_dim lat_dim time_dim]);
  stmp_var=netcdf.defVar(ncid2,'stmp','float',[lon_dim lat_dim time_dim]);
  spfh_var=netcdf.defVar(ncid2,'spfh','float',[lon_dim lat_dim time_dim]);
  netcdf.putAtt(ncid2,time_var,'base_date',int32([year mont day 0]));
  netcdf.endDef(ncid2);

  %rad
  lon_dim3=netcdf.defDim(ncid3,'nx_grid',nx);
  lat_dim3=netcdf.defDim(ncid3,'ny_grid',ny);
  time_dim3=netcdf.defDim(ncid3,'time',ntime);
  time_var3=netcdf.defVar(ncid3,'time','float',time_dim3);
  lon_var3=netcdf.defVar(ncid3,'lon','float',[lon_dim3 lat_dim3]);
  lat_var3=netcdf.defVar(ncid3,'lat','float',[lon_dim3 lat_dim3]);
  longwv_var=netcdf.defVar(ncid3,'dlwrf','float',[lon_dim3 lat_dim3 time_dim3]);
  solar_var=netcdf.defVar(ncid3,'dswrf','float',[lon_dim3 lat_dim3 time_dim3]);
  netcdf.putAtt(ncid3,time_var3,'base_date',int32([year mont day 0]));
  netcdf.endDef(ncid3);

  %prc
  lon_dim4=netcdf.defDim(ncid4,'nx_grid',nx);
  lat_dim4=netcdf.defDim(ncid4,'ny_grid',ny);
  time_dim4=netcdf.defDim(ncid4,'time',ntime);
  time_var4=netcdf.defVar(ncid4,'time','float',time_dim4);
  lon_var4=netcdf.defVar(ncid4,'lon','float',[lon_dim4 lat_dim4]);
  lat_var4=netcdf.defVar(ncid4,'lat','float',[lon_dim4 lat_dim4]);
  prc_var=netcdf.defVar(ncid4,'prate','float',[lon_dim4 lat_dim4 time_dim4]);
  netcdf.putAtt(ncid4,time_var4,'base_date',int32([year mont day 0]));
  netcdf.endDef(ncid4);

  %end define mode
  %write data - air
  netcdf.putVar(ncid2,lon_var,xgrid);
  netcdf.putVar(ncid2,lat_var,ygrid);
  netcdf.putVar(ncid2,time_var,single(time2/86400.));
  netcdf.putVar(ncid2,press_var,press);
  netcdf.putVar(ncid2,u10_var,u10);
  netcdf.putVar(ncid2,v10_var,v10);
  netcdf.putVar(ncid2,stmp_var,t2);%JP: add 2°C warming
  netcdf.putVar(ncid2,spfh_var,single(spfh));

  %Rad
  netcdf.putVar(ncid3,lon_var3,xgrid);
  netcdf.putVar(ncid3,lat_var3,ygrid);
  netcdf.putVar(ncid3,time_var3,single(time2/86400.));
  netcdf.putVar(ncid3,longwv_var,single(longwv)); %JP: add radiative forcing (warming)
  netcdf.putVar(ncid3,solar_var,single(solar));

  %prc
  netcdf.putVar(ncid4,lon_var4,xgrid);
  netcdf.putVar(ncid4,lat_var4,ygrid);
  netcdf.putVar(ncid4,time_var4,single(time2/86400.));
  netcdf.putVar(ncid4,prc_var,single(prc/dt)); %convert to rate

  if(1==2)
    for sl=1:ntime
      start=[0 0 sl-1]; %tcount];
      count=[nx ny 1];
      press2=squeeze(press(:,:,sl));
      netcdf.putVar(ncid2,press_var,start,count,press2);
      u10_2=squeeze(u10(:,:,sl));
      netcdf.putVar(ncid2,u10_var,start,count,u10_2);
      v10_2=squeeze(v10(:,:,sl));
      netcdf.putVar(ncid2,v10_var,start,count,v10_2);
      t2_2=squeeze(t2(:,:,sl));
      netcdf.putVar(ncid2,stmp_var,start,count,t2_2);
      netcdf.putVar(ncid2,spfh_var,start,count,single(spfh));
      netcdf.putVar(ncid2,time_var,sl-1,1,time2(sl)); %single(time2/86400.));
    end %sl
  end %if 1==2

  netcdf.close(ncid2);
  netcdf.close(ncid3);
  netcdf.close(ncid4);
end %d
return
% % % % %readnc2.m %Read in our sflux netcdf files into matlab and plot; modify as appro.
% % % % % e.g., the start and end stack #; time stamp for plot; sub-sampling freq.
% % % % % in plot
% % % % clear all; close all;
% % % % %scrsz = get(0,'ScreenSize'); %screen size
% % % % %Dimension for arrays reversed from ncdump (FORTRAN convention)
% % % % 
% % % % CB_bnd=load('CB_bnd.xy'); %load domain bnd
% % % % 
% % % % %NARR files
% % % % fill_in=1.e9; %junk value from nc files
% % % % avi_out = avifile('sflux.avi');
% % % % for i=1:10 %stack # for nc files
% % % %   char=sprintf('%3.3d',i);
% % % %   %filen=strcat('sflux_air_1.',char,'.nc');
% % % %   ncid0 = netcdf.open(filen,'NC_NOWRITE');
% % % %   vid=netcdf.inqVarID(ncid0,'time'); %input var./array name
% % % %   time_narr= netcdf.getVar(ncid0, vid);
% % % %   vid=netcdf.inqVarID(ncid0,'lon');
% % % %   lon_narr = netcdf.getVar(ncid0, vid); 
% % % %   vid=netcdf.inqVarID(ncid0,'lat');
% % % %   lat_narr = netcdf.getVar(ncid0, vid); 
% % % % 
% % % %   vid=netcdf.inqVarID(ncid0,'uwind');
% % % %   %Time dimension is last index in uwind_narr etc.
% % % %   uwind_narr= netcdf.getVar(ncid0, vid); 
% % % %   uwind_narr(find(abs(uwind_narr)>fill_in))=nan;
% % % %   vid=netcdf.inqVarID(ncid0,'vwind');
% % % %   vwind_narr= netcdf.getVar(ncid0, vid); 
% % % %   vwind_narr(find(abs(vwind_narr)>fill_in))=nan;
% % % %   vid=netcdf.inqVarID(ncid0,'prmsl');
% % % %   pres_narr= netcdf.getVar(ncid0, vid); 
% % % %   pres_narr(find(abs(pres_narr)>fill_in))=nan;
% % % %   vid=netcdf.inqVarID(ncid0,'stmp');
% % % %   airt_narr= netcdf.getVar(ncid0, vid); %ait temp.
% % % %   airt_narr(find(abs(airt_narr)>fill_in))=nan;
% % % %   vid=netcdf.inqVarID(ncid0,'spfh'); %humidity
% % % %   spfh_narr= netcdf.getVar(ncid0, vid); %specific humidity
% % % %   spfh_narr(find(abs(spfh_narr)>fill_in))=nan;
% % % % 
% % % %   disp(strcat('Done reading: ',filen));
% % % % 
% % % %   for j=1:length(time_narr)
% % % %     %Time stamp for plot - modify as appropr.
% % % %     date_narr=strcat('Sept ',num2str(i),': ',num2str(time_narr(j)*24),' hour UTC');
% % % %     %Compute bound
% % % %     xmin=min(min(lon_narr)); 
% % % %     ymin=min(min(lat_narr));
% % % %     xmax=max(max(lon_narr));
% % % %     ymax=max(max(lat_narr));
% % % % 
% % % %     %plot domain bnd
% % % %     hold on;
% % % %     plot(CB_bnd(:,1),CB_bnd(:,2),'k'); 
% % % % 
% % % %     %Scalar plot
% % % % %    contour(lon_narr,lat_narr,pres_narr(:,:,j));
% % % % %    caxis([9.7e4 10.2e4]);
% % % % %    colorbar;
% % % % 
% % % %     %Subsample
% % % %     ind1=1:1:size(lon_narr,1);
% % % %     ind2=1:1:size(lon_narr,2);
% % % % 
% % % %     %Vel. plot
% % % %     scale=0.2; %scale to fit
% % % %     quiver(lon_narr(ind1,ind2),lat_narr(ind1,ind2),...
% % % %     uwind_narr(ind1,ind2,j)*scale,vwind_narr(ind1,ind2,j)*scale,0,'b');
% % % %     quiver(0,56,5*scale,0.,0,'r');
% % % %     text(0,56.5,'5 m/s');
% % % % 
% % % %     %Atmos. pressure
% % % %     pcolor(lon_narr(ind1,ind2),lat_narr(ind1,ind2),double(pres_narr(ind1,ind2,j)/1.e5));
% % % %     shading interp; %to get rid of gridline
% % % %     colorbar;
% % % % 
% % % %     text(-2,56,date_narr);
% % % %     %Warning: avi does not like changes in xmin etc; use fixed scale
% % % %     axis([-10 15 30 60]);
% % % %     xlabel('Lon'); ylabel('Lat');
% % % %     title('Atmos. pressure scaled by 1.e5');
% % % % 
% % % % %   Stop here for testing
% % % % %    return;
% % % % 
% % % %     frame = getframe(gca);
% % % %     avi_out=addframe(avi_out,frame);
% % % %     clf; %clear figure
% % % %   end %j
% % % % 
% % % %   clear time_narr lon_narr lat_narr uwind_narr vwind_narr pres_narr airt_narr spfh_narr;
% % % % end %for all nc files
% % % % avi_out=close(avi_out);
%short_wave_radiation.m %  function [solar]=short_wave_radiation(yday,hour,lat,lon,tcc)
%  Calculate short wave (solar) radiation.
%  Adapted from Karsten Bolding and Hans Burchard
%  Short wave radiation is calculated based on the following input
%  parameters - year day, hour of day, latitude and longitude and cloud cover.
%  The albedo monthly values are from Payne (1972) as means of the values
%  at 40N and 30N for the Atlantic Ocean ( hence the same latitudinal
%  band of the Mediterranean Sea ) :
%  The radiation is returned as $W/m^2$.
%  Joseph's notes: should be valid gobally.
%  Inputs
%    integer  :: yday
%    double:: hour (UTC)
%    double arrays :: lat,lon,tcc (1D arrays) - tcc is cloud coverage (0-1)
%  Outputs: double array :: solar

  function [solar]=short_wave_radiation(yday,hour,lon,lat,tcc,cloud_shadow_eff)

%Check array dimensions
  if(length(lon) ~= length(lat) || length(lon) ~= length(tcc))
    error('short_wave_radiation: Arrays size mismatch');
  end
  ndim=length(lon);

  deg2rad=pi/180.; rad2deg=180./pi;
  solar=1350.;
  eclips=23.439*deg2rad;
  tau=0.7;
  aozone=0.09;
  yrdays=[365. 366.];
  alb1=[.719 .656 .603 .480 .385 .300 .250 .193 .164  .131 ...
        .103 .084 .071 .061 .054 .039 .036 .032 .031 .030];
  za=[90. 88. 86. 84. 82. 80. 78. 76. 74. 70.  ...
      66. 62. 58. 54. 50. 40. 30. 20. 10. 0.0];
  dza(1:8)=2.; dza(9:14)=4.; dza(15:19)=10.;

  th0 = 2*pi*yday/yrdays(1);
  th02 = 2*th0;
  th03 = 3*th0;

% sun declination :
  sundec = 0.006918 - 0.399912*cos(th0) + 0.070257*sin(th0) ... 
           - 0.006758*cos(th02) + 0.000907*sin(th02)       ... 
           - 0.002697*cos(th03) + 0.001480*sin(th03);

   alon = deg2rad*lon;
   alat = deg2rad*lat;

%  sun hour angle :
   thsun=(hour-12)*15*deg2rad+alon;

%  cosine of the solar zenith angle :
   coszen =sin(alat)*sin(sundec)+cos(alat)*cos(sundec).*cos(thsun);
   %Conserve memory
   clear thsun;
   [i1]=find(coszen>0);
   [i2]=find(coszen<=0);
   %To avoid overflow
   [i3]=find(coszen>0 & coszen<0.01);
   coszen(i3)=0.01;

   qatten(i2)=0;
   coszen(i2)=0;
  
   %contourf(coszen); 
   %colorbar;  

   qatten(i1)=tau.^(1./coszen(i1));

%   if(coszen <= 0.0) 
%     coszen = 0.0
%     qatten = 0.0
%   else
%     qatten = tau^(1./coszen);
%   end 
   qzer=coszen*solar;
   qdir=qzer.*qatten'; %qatten is transposed
   clear qatten;
   qdiff=((1-aozone)*qzer-qdir)*0.5;
   clear qzer;
   qtot=qdir+qdiff;
   clear qdiff qdir;
   eqnx=(yday-81)/yrdays(1)*2*pi;

%  sin of the solar noon altitude in radians 
   sunbet=sin(alat)*sin(eclips*sin(eqnx))+cos(alat)*cos(eclips*sin(eqnx));
%  solar noon altitude in degrees 
   sunbet = asin(sunbet)*rad2deg;

%  calculates the albedo as a function of the solar zenith angle
%  (after Payne jas 1972)
%  solar zenith angle in degrees
   zen=(180./pi)*acos(coszen);
   [i1]=find(zen>=74);
   [i2]=find(zen>=50 & zen<74);
   [i3]=find(zen<50);
   jab(i1)=ceil(.5*(90.-zen(i1))+1);
   jab(i2)=ceil(.23*(74.-zen(i2))+9);
   jab(i3)=ceil(.10*(50.-zen(i3))+15);

%   if(zen .ge. 74.)then
%      jab=.5*(90.-zen)+1.
%   else if (zen .ge. 50.) then
%      jab=.23*(74.-zen)+9.
%   else
%      jab=.10*(50.-zen)+15.
%   endif

   indx=find(jab>19);
   if(length(indx) ~=0); error('short_wave_radiation: index out of bound'); end;

   %Have to use serial mode here for differencing
   for i=1:ndim
     dzen=(za(jab(i))-zen(i))/dza(jab(i));
     albedo(i)=alb1(jab(i))+dzen*(alb1(jab(i)+1)-alb1(jab(i)));
   end %for i
   rat=(1-cloud_shadow_eff*tcc+.0019*sunbet).*(1-albedo'); %rat=(1-0.62*tcc+.0019*sunbet).*(1-albedo'); original equaion
   [i1]=find(rat>1);
   rat(i1)=1;
   solar=qtot.*rat;

return
