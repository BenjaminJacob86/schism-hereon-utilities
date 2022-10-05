function gen_ncBDforcing
%
% generates netcdf boundary forcing for SCHISM
% based on different amm7 products
%
% 1) perform horizontally inversed distance (with respect to ocean nodes)
%    interpolation to get profiles at SCHISM boundary nodes
% 2) perform linear vertical interpolation on those profiles to yield
%    profiles at SCHISM node depth
%  
%    the sea level fluctations are not accounted for in computation
%    of SCHISM zcoordinates i.e. tracer are interpolated assuming
%    zero sea level. but error is likely small with respect
%    to vertical resolution of data anyways.
%
%IN: netcdf, sowie hgrid.ll mit den Punkten und Randpolygonen, vgrid.in
%OUT: binary *.th
%
% author: Benjamnin Jacob 01.07.2019
%
% 
% 22.06.2020  change  hardcoded date suffix to * and dir seach for ile
% since distance between hindcast day and simulated day may very in names of amm input

close all;
clc;

%% Parameter - - Userinterface
%
% Simulation StartDate and Period
%startdate-enddate: target period
startdate=datenum(2017,01,02);  %my ocean begins with 1h and ends with 00
restartdate=startdate; %datenum(2016,03,02) % zero padding at beginning if restartdate != startdate; if forcing output is intended to continue after restart.
enddate=datenum(2018,12,31);
ndays=enddate-restartdate+1; % number of days

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Directories %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rundir='/gpfs/work/jacobb/data/RUNS/GermanBight/GermanBight_2017_2018/'; %The Model Setup Directory
bnddir=rundir; % Directory where the boundaryforcing is to be saved
sourcedir='/gpfs/work/ksddata/model/ocean/AMM15/';
addpath(rundir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

donedays=restartdate-startdate;
timestepadd=donedays*86400;

%%%%%%%%%%%%%%% Inputfiles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% choose data set  via amm7_format
amm7_format=4;  % 1: my ocean 3 layers original names 2: my ocean 3 layers  3: my Ocean 24 layers  4: Amm 15

sufix_elev_filename=('.nc');  %elevation
sufix_salt_filename=('.nc');  %salt
sufix_temp_filename=('.nc');  %temperature
sufix_vel_filename=('.nc');   %velocity


recheck_varids=1;              % recheck varids each time a file is load ( I experienced chaning ids within 2018 for copernicus data)

switch amm7_format
    case 1 % Use My Ocean 3 original names
        vrt_in_sigma=0;  % vertical coordinate is (0) Z (1) sigma -> different interpol
        downard_sign=1; % only sigma coordinates; SELFE negative downward
        
        prefix_elev_filename=('metoffice_foam1_amm7_NWS_ASLV'); %elevation
        prefix_salt_filename=('metoffice_foam1_amm7_NWS_SAL');  %salt
        prefix_temp_filename=('metoffice_foam1_amm7_NWS_TEMP'); %temperature
        prefix_vel_filename=('metoffice_foam1_amm7_NWS_CUR');   %velocity
        
        %datetag generation function
        filename_datetag=@(day_date)[ '_b' datestr(day_date+1,'yyyymmdd')  '_hi' datestr(day_date,'yyyymmdd')];
        
        label='Amm7_3old';
        FillValue=-32768;
    case 2  % Use My Ocean 3 layer
        
        vrt_in_sigma=1;  % vertical coordinate is (0) Z (1) sigma -> different interp
        downard_sign=1; % only sigma coordinates; SELFE negative downward
        % sigmalayers multiplied by -1 * downward sign
        prefix_elev_filename=('metoffice_foam1_amm7_NWS_ASLV');     %elevation
        prefix_salt_filename=('metoffice_foam1_amm7_NWS_PSAL');  %salt
        prefix_temp_filename=('metoffice_foam1_amm7_NWS_TEMP');  %temperature
        prefix_vel_filename=('metoffice_foam1_amm7_NWS_RFVL');   %velocity
        
        %datetag generation function
        %filename_datetag=@(day_date)[datestr(day_date,'yyyymmdd')];
        filename_datetag=@(day_date)[ '_b' datestr(day_date+1,'yyyymmdd')  '_hi' datestr(day_date,'yyyymmdd')];
        
        label='Amm7_3';
        FillValue=-32768;
    case 3 % UseA,,7  My Ocean 24
        
        vrt_in_sigma=0;  % vertical coordinate is (0) Z (1) sigma -> different interpol
        downard_sign=1; % only sigma coordinates; SELFE negative downward
        % sigmalayers multiplied by -1 * downward sign
        %prefix_elev_filename=('metoffice_foam1_amm7_NWS_SSH');     %elevation
        %prefix_salt_filename=('metoffice_foam1_amm7_NWS_SAL');  %salt
        %prefix_temp_filename=('metoffice_foam1_amm7_NWS_TEM');  %temperature
        %prefix_vel_filename=('metoffice_foam1_amm7_NWS_CUR');   %velocity   %velocity
        prefix_elev_filename=('metoffice_foam1_amm7_NWS_ASLV');     %elevation
        prefix_salt_filename=('metoffice_foam1_amm7_NWS-HR_PSAL');  %salt
        prefix_temp_filename=('metoffice_foam1_amm7_NWS-HR_TEMP');  %temperature
        prefix_vel_filename=('metoffice_foam1_amm7_NWS-HR_RFVL');   %velocity
        
        %datetag generation function
        filename_datetag=@(day_date)[ '_b' datestr(day_date+1,'yyyymmdd')  '_hi' datestr(day_date,'yyyymmdd')];
        label='Amm7_24';
        FillValue=-32768;
        
        
    case 4 % UseA,,7  My Ocean 24
        
        vrt_in_sigma=0;  % vertical coordinate is (0) Z (1) sigma -> different interpol
        downard_sign=1; % only sigma coordinates; SELFE negative downward
        % sigmalayers multiplied by -1 * downward sign
        prefix_elev_filename=('metoffice_foam1_amm15_NWS_SSH');     %elevation
        prefix_salt_filename=('metoffice_foam1_amm15_NWS_SAL');  %salt
        prefix_temp_filename=('metoffice_foam1_amm15_NWS_TEM');  %temperature
        prefix_vel_filename=('metoffice_foam1_amm15_NWS_CUR');   %velocity   %velocity
        
        %datetag generation function
		% not necesary fixey
        filename_datetag=@(day_date)[ '_b' datestr(day_date+1,'yyyymmdd')  '_hi' datestr(day_date,'yyyymmdd')];
		filename_datetag=@(day_date)[ '_b*_hi' datestr(day_date,'yyyymmdd')];
        label='Amm15';
        FillValue=-32768;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Outputfiles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read vertical layers from vgrid in
fid=fopen(fullfile(rundir, 'vgrid.in' ),'r');
line=fgets(fid);
ivcor=textscan(line,'%f',1);ivcor=ivcor{1};
if ivcor==1 % LSC 2
    flexible_vertical_coordinate=1;
elseif ivcor==2
    flexible_vertical_coordinate=0;
end
line=fgets(fid);
nrvt=textscan(line,'%f',1);nrvt=nrvt{1};
fclose(fid);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Interp Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose the method you want to use to create the boundary conditions:
% 'old_method' original method, 'lin_interp' linear interpolation, 'next_neighbour' next neighbour

vert_interp_meth='linear';%'linear','nearest'
hor_interp_meth=1; %0: 'nearest', 1: 'linear' ; % so far only nearest neighbour implementet
if hor_interp_meth==0
    horintstr='nearest'
else
    horintstr='linterp'
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

suffix='_nlslayer_'; %namens anhang
if startdate==restartdate
    daterange=['_' datestr(startdate,'YYmmdd') '-'  datestr(enddate,'YYmmdd')];
else
    daterange=['_' datestr(startdate,'YYmmdd') '_' datestr(restartdate,'YYmmdd') '-'  datestr(enddate,'YYmmdd')];
end
elevFileOut=['elev3D_' horintstr suffix daterange  '.th'];     %elevation
saltFileOut=['salt3D_' horintstr suffix daterange  '.th'];     %salt
tempFileOut=['temp3D_' horintstr suffix daterange  '.th'];     %temperature
uvFileOut=  ['uv3D_' horintstr suffix daterange    '.th'];       %u v component current
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



disp(['using forcing data '  label])
%% Bestimme open Boundary aus Boundaries in hgrid.gr3 und forcing spezifikation in bctides.in
disp('Bestimme Forcingr???nder')
disp(['Forcing is in sigma layers: ' num2str(vrt_in_sigma)])
disp(['horizontal grid interpolation: ' hor_interp_meth ])
disp(['Vertical profile interpolation: ' vert_interp_meth ])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Load information from hgrid.gr3 %%%%%%%%%%%%%%%%%%%%%%%%
fid=fopen(fullfile(rundir,'hgrid.gr3'),'r');
textline=fgets(fid);
pattern=' Number of open boundaries';
while ~contains(textline,pattern)
    textline=fgets(fid);
end

n_opbd=textscan(textline,'%d ');n_opbd=n_opbd{1};
fgets(fid);%textline=fgets(fid);
%n_opbd_nodes=textscan(textline,'%d ');n_opbd_nodes=n_opbd_nodes{1};
openbds.Nopenbd=n_opbd;

for i=1:n_opbd
    textline=fgets(fid);%Number of nodes for open bundary i
    temp=textscan(textline,'%d ');
    n_nodes=temp{1};
    openbds.boundaries{i}.n_nodes=n_nodes;
    openbds.boundaries{i}.nodes=textscan(fid,'%d \n',n_nodes);%read boundaries segments (nr of node specified by bound segs(i)
end
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Load information from bctides.in%%%%%%%%%%%%%%%%%%%%%%%%
searchfrom='%d %d %d %d %d'; %Suche Zeilen mit format Nboundaries 4 4 4 4
% 4 is 3D timehistory forcing
fid=fopen(fullfile(bnddir,'bctides.in'),'r');

line=fgets(fid);
k=0;   % counter for 3D TH forced Boundaries
while ~feof(fid)
    linevals=textscan(line,searchfrom);
    %Bctides Zeile erf???llt format  und hat bei einer der variablen eine 4 als
    %mareker f???r time history forcing
    %JP: eigentlich genuegt wenn Salz oder T eine 4 haben
    if ~isempty(linevals{5}) && (linevals{4}==4  && linevals{5}==4 )% || linevals{4}==4  || linevals{5}==4 )
        k=k+1;
        %ans=linevals{1}; %JP
        %THboundary.nnodes{k}=ans(1);
        THboundary.nnodes{k}=linevals{1}(1);
    end
    line=fgets(fid);
end
fclose(fid);


bdlist=1:openbds.Nopenbd;

i_thbounds=[];  %Nummern der offenen r???ndern an Denenen Timehistory forcing gesetzt ist
openbd=[];

Len_force_bd=[];
for thbound=1:k
    nmember=0;
    for i=bdlist
        openbds.boundaries{i}.n_nodes;
        % Nr boundary nodes von boundarie f???r timehistory forcing entspricht der
        % Nummer des aktullen randes ? -> rand gefunden
        nmember=nmember+1;
        if openbds.boundaries{i}.n_nodes==THboundary.nnodes{thbound}
            %stop
            bdlist(nmember)=[];
            i_thbounds=[i_thbounds i];
            openbd=[openbd; openbds.boundaries{i}.nodes{1}];
            Len_force_bd=[Len_force_bd openbds.boundaries{i}.n_nodes];
            break
        end
    end
end
fprintf('Forcing at %i openboundaries:',numel(Len_force_bd))
fprintf('Open Boundary %i with %i nodes\n',[i_thbounds; Len_force_bd])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% knoten alle offenen grenzen f???r forcing sind nun in openbd vereint
%%  Lade Modell grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Lade Modellgitter von Zielmodell %%%%%%%% %%%%%%%%%%%%%%%%%
idG=fopen(fullfile(rundir,'hgrid.ll'));
fgets(idG);%Grid name
nnde=textscan(idG,'%d %d',1);%Anzahl Elemente | Anzahl nodes
nnde=nnde{2};
form='%f %f %f %f'; %Format des Node abschnitts
modelgrid=textscan(idG,form,nnde);
fclose(idG);
modelgrid=cell2mat(modelgrid);
bd_salt=modelgrid(openbd,:);
nbd=length(openbd); %length of open boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r=bd_salt(:,4)<0;
bd_salt(r,4)=10; %JP: f?r trockene Randknoten 10 m Tiefe annehmen


%% start interpolation routine

t_total=0;  %total timestep counter
n_total=1+ndays*24; % all days and one from step before
%% init i/o

% intput
dy=-1;
date_suffix=filename_datetag(restartdate+dy);
filessh=dir(fullfile(sourcedir,[prefix_elev_filename date_suffix sufix_elev_filename]));
filesalt=dir(fullfile(sourcedir,[prefix_salt_filename date_suffix sufix_salt_filename]));
filetemp=dir(fullfile(sourcedir,[prefix_temp_filename date_suffix sufix_temp_filename]));
filevel=dir(fullfile(sourcedir,[prefix_vel_filename date_suffix sufix_vel_filename]));

ncids=struct('elev',netcdf.open(fullfile(filessh.folder,filessh.name),'NC_NOWRITE'),...
    'salt',netcdf.open(fullfile(filesalt.folder,filesalt.name),'NC_NOWRITE'),...
    'temp',netcdf.open(fullfile(filetemp.folder,filetemp.name),'NC_NOWRITE'),...
    'vel',netcdf.open(fullfile(filevel.folder,filevel.name),'NC_NOWRITE'));


varids=struct('elev',...
    struct('data',get_varid_from_standard_name(ncids.elev,'sea_surface_height_above_geoid'),... % varids.elev=netcdf.inqVarID(ncids.elev,'sossheig');
    'time',netcdf.inqVarID(ncids.elev,'time')),...
    'salt',struct(...
    'data',get_varid_from_standard_name(ncids.salt,'sea_water_salinity'),... %     varids.salt=netcdf.inqVarID(ncids.salt,'vosaline');
    'depth',netcdf.inqVarID(ncids.salt,'depth')),...
    'temp',struct(...
    'data',get_varid_from_standard_name(ncids.temp,'sea_water_potential_temperature'),... %     varids.salt=netcdf.inqVarID(ncids.salt,'vosaline');
    'depth',netcdf.inqVarID(ncids.temp,'depth')),...
    'u',struct(...
    'data',get_varid_from_standard_name(ncids.vel,'eastward_sea_water_velocity'),... %     varids.salt=netcdf.inqVarID(ncids.salt,'vosaline');
    'depth',netcdf.inqVarID(ncids.salt,'depth')),...
    'v',struct(...
    'data',get_varid_from_standard_name(ncids.vel,'northward_sea_water_velocity'))... %     varids.salt=netcdf.inqVarID(ncids.salt,'vosaline');
    );

%


%%
% check time step of input data
time_elev=double(netcdf.getVar(ncids.elev,varids.elev.time));
unitstr=netcdf.getAtt(ncids.elev,varids.elev.time,'units');
if contains(unitstr,'seconds')
    str0=strfind(unitstr,'since ')+1+length('since');
    forcing_dt=diff(time_elev(1:2));
    selfe_dt=forcing_dt;
    reftimein=datenum(unitstr(str0:end));
else
    error('forcing data time not in sceonds Adjust code')
end
nt=ndays*86400/forcing_dt+1;
time_nc=zeros(nt,1);
ntin=length(time_elev);
t0=ntin;
timein=zeros(ntin*ndays+1,1);
offset_time=((reftimein+time_elev(t0)/86400)-startdate)*86400;  % offset between first used time step of input data and desired timesteps (used for interpolation)

%% get scale factors
%if strcmp(label,'Amm7_24') | strcmp(label,'Amm15') % read scale factort from netcdf
scale_factor_elev=double(netcdf.getAtt(ncids.elev,varids.elev.data,'scale_factor'));
scale_factor_salt=double(netcdf.getAtt(ncids.salt,varids.salt.data,'scale_factor'));
scale_factor_temp=double(netcdf.getAtt(ncids.temp,varids.temp.data,'scale_factor'));
scale_factor_vel=netcdf.getAtt(ncids.vel,varids.u.data,'scale_factor');

%offset_elev=double(netcdf.getAtt(ncids.elev,0,'add_offset'));
offset_temp=double(netcdf.getAtt(ncids.temp,varids.temp.data,'add_offset'));
offset_salt=double(netcdf.getAtt(ncids.salt,varids.salt.data,'add_offset'));
%offset_vel=netcdf.getAtt(ncids.vel,0,'add_offset');

% include conversion from kelvin to degree in offset_temp, if data is in kelvin
if strcmpi(netcdf.getAtt(ncids.temp,varids.temp.data,'units'),'k') || strcmpi(netcdf.getAtt(ncids.temp,varids.temp.data,'units'),'kelvin')
    offset_temp=offset_temp-273.15;
end

%% load data grid
varids.elev_lat=netcdf.inqVarID(ncids.elev,'lat');
varids.elev_lon=netcdf.inqVarID(ncids.elev,'lon');
lat_elev=double(netcdf.getVar(ncids.elev,varids.elev_lat));
lon_elev=double(netcdf.getVar(ncids.elev,varids.elev_lon));


varids.salt_lat=netcdf.inqVarID(ncids.salt,'lat');
varids.salt_lon=netcdf.inqVarID(ncids.salt,'lon');
lat_salt=double(netcdf.getVar(ncids.salt,varids.salt_lat));
lon_salt=double(netcdf.getVar(ncids.salt,varids.salt_lon));

varids.temp_lat=netcdf.inqVarID(ncids.temp,'lat');
varids.temp_lon=netcdf.inqVarID(ncids.temp,'lon');
lat_temp=double(netcdf.getVar(ncids.temp,varids.temp_lat));
lon_temp=double(netcdf.getVar(ncids.temp,varids.temp_lon));

varids.vel_lat=netcdf.inqVarID(ncids.vel,'lat');
varids.vel_lon=netcdf.inqVarID(ncids.vel,'lon');
lat_vel=double(netcdf.getVar(ncids.vel,varids.vel_lat));
lon_vel=double(netcdf.getVar(ncids.vel,varids.vel_lon));

[grids.elev.lat,grids.elev.lon] = meshgrid(lat_elev,lon_elev);
[grids.salt.lat,grids.salt.lon] = meshgrid(lat_salt,lon_salt);
[grids.temp.lat,grids.temp.lon] = meshgrid(lat_temp,lon_temp);
[grids.vel.lat,grids.vel.lon] = meshgrid(lat_vel,lon_vel);


%% determin data dimensios and necessary grid range
%
elev_in=double(netcdf.getVar(ncids.elev,varids.elev.data));
time_elev=double(netcdf.getVar(ncids.elev,varids.elev.time));
salt_in=double(netcdf.getVar(ncids.salt,varids.salt.data)); % Fillvalue FillValue
[a, b, nz, ntime]=size(salt_in);
temp_in=double(netcdf.getVar(ncids.temp,varids.temp.data)); % Fillvalue FillValue
u_in=double(netcdf.getVar(ncids.vel,varids.u.data)); % Fillvalue FillValue (lon,lat,depth,time)
v_in=double(netcdf.getVar(ncids.vel,varids.v.data)); % Fillvalue FillValue (lon,lat,depth,time)
dims={size(elev_in),size(salt_in),size(temp_in),size(u_in)};
%clear elev_in, salt_in, temp_in, u_in

%% limit grid to later only load subset of ncfile
%% get sorounding nodes on all variable grids
% parents
parentsii=zeros(length(bd_salt),2,5); % 4 structured grid nodes sorounding destination poibts
parentsjj=parentsii; % for all 4 variables [openbd node,index,variable]
% 1: elev 2:salt 3:temp 4:u 5: v
lons={lon_elev,lon_salt,lon_temp,lon_vel,lon_vel};
lats={lat_elev,lat_salt,lat_temp,lat_vel,lat_vel};
%dists=zeros(2,2,length(bd_salt),5);
ws=zeros(2,2,length(bd_salt),nz,5); %inverse distance weights

for i=1:size(bd_salt,1) %JP: 1->2
    % read the target position
    lon_p = bd_salt(i,2);
    lat_p = bd_salt(i,3);
    for ivar=1:5
        % find the four closest nodes from rect grid
        [~,ii]=sort(abs(lons{ivar}-lon_p));ii=sort(ii(1:2)); % first dim
        [~,jj]=sort(abs(lats{ivar}-lat_p));jj=sort(jj(1:2)); % second dim
        
        % -> lat
        % |
        % v lon
        % bottom left   top left
        % bottom right  top right
        %size(elev)==[length(lon_elev) length(lat_elev)]
        % [YY,XX]=meshgrid(lat_elev,lon_elev); % has to be same dimension as variable matrox elev etc.
        %  XX(ii,jj), YY(ii,jj)
        % before bilinear interpolation replace parant indieces of dry nodes with closest non dry node indices
        switch ivar
            case 1
                ifill=elev_in(ii,jj)==FillValue;
                loni=grids.elev.lon(ii,jj);
                lati=grids.elev.lat(ii,jj);
            case 2
                ifill=salt_in(ii,jj,:,1)==FillValue;
                loni=grids.salt.lon(ii,jj);
                lati=grids.salt.lat(ii,jj);
            case 3
                ifill=temp_in(ii,jj,:,1)==FillValue;
                loni=grids.temp.lon(ii,jj);
                lati=grids.temp.lat(ii,jj);
            case 4
                ifill=u_in(ii,jj,:,1)==FillValue;
                loni=grids.vel.lon(ii,jj);
                lati=grids.vel.lat(ii,jj);
            case 6
                ifill=v_in(ii,jj,:,1)==FillValue;
                loni=grids.vel.lon(ii,jj);
                lati=grids.vel.lat(ii,jj);
        end
        
        
        %compute prjection distance in cartesion coordinades
        [x,y]=ll2xy(lati,loni,lat_p,lon_p); % from batry
        %D=sqrt((loni-lon_p).^2+(lati-lat_p).^2); lon/lat distances
        D=sqrt(x.^2+y.^2); %
        
        % use matlas 1/inf = 0  to ignore weight of dry points
        %  thhis might be depth dependent also
        if ivar==1
            D(ifill)=inf;
            D=1./D;
            ws(:,:,i,1,ivar)=D./sum(D(:));%1./D/sum(1./D(:)); % take dry valuaes out of weights
        else
            D=repmat(D,1,1,nz);
            D(ifill)=inf;
            D=1./D;
            dsum=sum(sum(D,1),2);
            ws(:,:,i,:,ivar)=D./dsum;%1./D/sum(1./D(:)); % take dry valuaes out of weights
        end
        
        %#
        %make inverse distance weight of non dry eleents
        
        parentsii(i,:,ivar)=ii;
        parentsjj(i,:,ivar)=jj;
    end
end


% check at end of boundaries if inlet and move accoring le
% double at bounary
% expand left boundaries
i1=1;
for ibnd = 1:length(THboundary.nnodes)
        i2=i1+THboundary.nnodes{ibnd}-1;
        if max(max(max(max(ws(:,:,i1,:,:)))))==0
            parentsii(i1,:,:)=parentsii(i1+1,:,:);
            parentsjj(i1,:,:)=parentsjj(i1+1,:,:);
            ws(:,:,i1,:,:)=ws(:,:,i1+1,:,:);
            ws(:,:,i1,:,:)=ws(:,:,i1+1,:,:);
            
        end
        if max(max(max(max(ws(:,:,i2,:,:)))))==0
            parentsii(i2,:,:)=parentsii(i2-1,:,:);
            parentsjj(i2,:,:)=parentsjj(i2-1,:,:);
            ws(:,:,i2,:,:)=ws(:,:,i2+1,:,:);
            ws(:,:,i2,:,:)=ws(:,:,i2+1,:,:);
        end
        i1=i1+THboundary.nnodes{ibnd};
        
end


% index range of netcdf covering domain
% to only read area insidy boundary boundy box from netcdf
iirange=[min(parentsii(:)) max(parentsii(:))];
jjrange=[min(parentsjj(:)) max(parentsjj(:))];

for ivar=1:4
    start{ivar}=zeros(1,length(dims{ivar}));
    count{ivar}=dims{ivar};
    start{ivar}(dims{ivar}==length(lons{ivar}))=iirange(1)-1;
    start{ivar}(dims{ivar}==length(lats{ivar}))=jjrange(1)-1;
    count{ivar}(dims{ivar}==length(lons{ivar}))=diff(iirange)+1;
    count{ivar}(dims{ivar}==length(lats{ivar}))=diff(jjrange)+1;
end


grids.elev.lat=grids.elev.lat(iirange(1):iirange(2),jjrange(1):jjrange(2));
grids.elev.lon=grids.elev.lon(iirange(1):iirange(2),jjrange(1):jjrange(2));
grids.salt.lat=grids.salt.lat(iirange(1):iirange(2),jjrange(1):jjrange(2));
grids.salt.lon=grids.salt.lon(iirange(1):iirange(2),jjrange(1):jjrange(2));
grids.temp.lat=grids.temp.lat(iirange(1):iirange(2),jjrange(1):jjrange(2));
grids.temp.lon=grids.temp.lon(iirange(1):iirange(2),jjrange(1):jjrange(2));
grids.vel.lat=grids.vel.lat(iirange(1):iirange(2),jjrange(1):jjrange(2));
grids.vel.lon=grids.vel.lon(iirange(1):iirange(2),jjrange(1):jjrange(2));

parentsii=parentsii-iirange(1)+1;
parentsjj=parentsjj-jjrange(1)+1; % range of really needed forcing grid
% plot(bd_salt(:,2),bd_salt(:,3));
% hold on
% for i=1:length(bd_salt(:,2))
% plot(grids.elev.lon(parentsii(i,:,1),parentsjj(i,:,1)),grids.elev.lat(parentsii(i,:,1),parentsjj(i,:,1)),'r+')
% end

depth_salt=double(netcdf.getVar(ncids.salt,varids.salt.depth)); %myocean hires nls z levels, ie interpolate from nlsz levels to 31 sigma layers
depth_temp=double(netcdf.getVar(ncids.temp,varids.temp.depth)); %myocean hires nls z levels, ie interpolate from nlsz levels to 31 sigma layers
depth_vel=double(netcdf.getVar(ncids.vel,varids.u.depth)); %myocean hires nls z levels, ie interpolate from nlsz levels to 31 sigma layers

%% preallocate output arrays for n data steps
% before temporal ibterpolation
valsarray_elev_nc=zeros(1,1,nbd,n_total);
valsarray_salt_nc=zeros(1,nrvt,nbd,n_total);
valsarray_temp_nc=valsarray_salt_nc;
valsarray_vel_nc=zeros(2,nrvt,nbd,n_total);


% inout data z frame - uniform in grided data
z_salt=depth_salt*-downard_sign;
z_temp=depth_temp*-downard_sign;
z_vel=depth_vel*-downard_sign;


% load schism depth levels at nodes
zsselfe=nan(nrvt,nbd);
if flexible_vertical_coordinate==0
    for i=1:size(bd_salt,1)
        if vrt_in_sigma==0  % uniform depth levels
            [sigma_selfe]=calculate_z_coordinate_selfe(bd_salt(i,4),nrvt);
            zsselfe(:,i)=bd_salt(i,end)*sigma_selfe;
        else   %b) % Vertical coordinate is already sigma
            zsselfe(:,i)=calculate_z_coordinate_selfe(bd_salt(i,4),nrvt);
        end
    end
    
    
else % new vertical coordinate
    ibottom=nan(1,nbd);
    fid=fopen(fullfile(rundir, 'vgrid.in' ),'r');
    line=fgets(fid);
    ivcor=textscan(line,'%f',1);ivcor=ivcor{1};
    if ivcor==1 % LSC 2
        flexible_vertical_coordinate=1;
    elseif ivcor==2
        flexible_vertical_coordinate=0;
    end
    line=fgets(fid);
    nrvt=textscan(line,'%f',1);nrvt=nrvt{1};
    
    
    nde=0;
    icount=0;
    while ~feof(fid) & nde <= max(openbd)
        nde=nde+1;
        line=fgets(fid);
        if ismember(nde,openbd)
            icount=icount+1;
            nrs=cell2mat(textscan(line,'%f'));
            zsselfe(nrs(2):nrvt,icount)=nrs(3:end);
            %zsselfe(1:nrs(2)-1,count)=-1;
            ibottom(icount)=nrs(2);
        end
    end
    fclose(fid);
    for i=1:size(bd_salt,1)
        if vrt_in_sigma==0  % uniform depth levels
            zsselfe(:,i)=bd_salt(i,end)*zsselfe(:,i); %sigma_selfe;
        end
    end
end


for dy=-1:ndays-1 % Loop over Timesteps
    %fprintf('prepare forcing for day %i / %i',dy+1,ndays)
    disp(['prepare forcing for day ',num2str(dy+1),' / ',num2str(ndays)])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% Load Netcdf data %%%%%%%%%%%%
    %date_suffix=filename_datetag(restartdate+dy);
    %ncids.elev=netcdf.open(fullfile(sourcedir,[prefix_elev_filename date_suffix %sufix_elev_filename]),'NC_NOWRITE');
    % ncids.salt=netcdf.open(fullfile(sourcedir,[prefix_salt_filename date_suffix sufix_salt_filename]),'NC_NOWRITE');
    %ncids.temp=netcdf.open(fullfile(sourcedir,[prefix_temp_filename date_suffix %sufix_temp_filename]),'NC_NOWRITE');
    %ncids.vel=netcdf.open(fullfile(sourcedir,[prefix_vel_filename date_suffix sufix_vel_filename]),'NC_NOWRITE');

	date_suffix=filename_datetag(restartdate+dy);
	filessh=dir(fullfile(sourcedir,[prefix_elev_filename date_suffix sufix_elev_filename]));
	filesalt=dir(fullfile(sourcedir,[prefix_salt_filename date_suffix sufix_salt_filename]));
	filetemp=dir(fullfile(sourcedir,[prefix_temp_filename date_suffix sufix_temp_filename]));
	filevel=dir(fullfile(sourcedir,[prefix_vel_filename date_suffix sufix_vel_filename]));

	ncids=struct('elev',netcdf.open(fullfile(filessh.folder,filessh.name),'NC_NOWRITE'),...
		'salt',netcdf.open(fullfile(filesalt.folder,filesalt.name),'NC_NOWRITE'),...
		'temp',netcdf.open(fullfile(filetemp.folder,filetemp.name),'NC_NOWRITE'),...
		'vel',netcdf.open(fullfile(filevel.folder,filevel.name),'NC_NOWRITE'));    
    
    if dy==-1 % load temporally invariant information
    else
        t0=1;
    end
    
    % change in ncids?
    if recheck_varids
        varids=struct('elev',...
    struct('data',get_varid_from_standard_name(ncids.elev,'sea_surface_height_above_geoid'),... % varids.elev=netcdf.inqVarID(ncids.elev,'sossheig');
    'time',netcdf.inqVarID(ncids.elev,'time')),...
    'salt',struct(...
    'data',get_varid_from_standard_name(ncids.salt,'sea_water_salinity'),... %     varids.salt=netcdf.inqVarID(ncids.salt,'vosaline');
    'depth',netcdf.inqVarID(ncids.salt,'depth')),...
    'temp',struct(...
    'data',get_varid_from_standard_name(ncids.temp,'sea_water_potential_temperature'),... %     varids.salt=netcdf.inqVarID(ncids.salt,'vosaline');
    'depth',netcdf.inqVarID(ncids.temp,'depth')),...
    'u',struct(...
    'data',get_varid_from_standard_name(ncids.vel,'eastward_sea_water_velocity'),... %     varids.salt=netcdf.inqVarID(ncids.salt,'vosaline');
    'depth',netcdf.inqVarID(ncids.salt,'depth')),...
    'v',struct(...
    'data',get_varid_from_standard_name(ncids.vel,'northward_sea_water_velocity'))... %     varids.salt=netcdf.inqVarID(ncids.salt,'vosaline');
    );
    end
    %
    % load data - only in ncesseary reatngle
    elev_in=double(netcdf.getVar(ncids.elev,varids.elev.data,start{1},count{1}));
    time_elev=double(netcdf.getVar(ncids.elev,varids.elev.time));
    salt_in=double(netcdf.getVar(ncids.salt,varids.salt.data,start{2},count{2})); % Fillvalue FillValue
    %time_salt=double(netcdf.getVar(ncids.salt,varids.salt_time));
    temp_in=double(netcdf.getVar(ncids.temp,varids.temp.data,start{3},count{3})); % Fillvalue FillValue
    %time_temp=double(netcdf.getVar(ncids.temp,varids.temp_time));
    u_in=double(netcdf.getVar(ncids.vel,varids.u.data,start{4},count{4})); % Fillvalue FillValue (lon,lat,depth,time)
    v_in=double(netcdf.getVar(ncids.vel,varids.v.data,start{4},count{4})); % Fillvalue FillValue
    %time_vel=time_salt ;%double(netcdf.getVar(ncids.vel,varids.vel_time)); JP
    
    
    % get for parent nodes for each grid cell
    
    
    for t=t0:ntime %time loop within file
        t_total=t_total+1;
        
        %timein(t_total)=time_elev(t);
        
        %interp in horizontal plane inverse distance weighted
        % land nodes of structured grid are taken out of weighting
        % assumption is position of Fillvalues done change horitontally
        % over time within forcing data grid
        
        for i=1:size(bd_salt,1) % bd node loop
            
            % interp horizontally inverse distance weighted along horizontal, to get profiles at bd_nodes
            
            elev=elev_in(parentsii(i,:,1),parentsjj(i,:,1),t).* ws(:,:,i,1,1)*scale_factor_elev;  % dummys already weighetd zero
            elev_interp=sum(elev(:));
            
            % pfofiles from input data
            salt=(salt_in(parentsii(i,:,1),parentsjj(i,:,1),:,t)*scale_factor_salt+offset_salt) .* squeeze(ws(:,:,i,:,2));  % dummys already weighetd zero
            salt_profile=squeeze(sum(sum(salt)));
            
            temp=(temp_in(parentsii(i,:,1),parentsjj(i,:,1),:,t)*scale_factor_temp+offset_temp) .* squeeze(ws(:,:,i,:,3));  % dummys already weighetd zero
            temp_profile=squeeze(sum(sum(temp)));
            
            u=(u_in(parentsii(i,:,1),parentsjj(i,:,1),:,t)*scale_factor_vel) .* squeeze(ws(:,:,i,:,4));  % dummys already weighetd zero
            u_profile=squeeze(sum(sum(u)));
            
            v=(v_in(parentsii(i,:,1),parentsjj(i,:,1),:,t)*scale_factor_vel) .* squeeze(ws(:,:,i,:,4));  % dummys already weighetd zero
            v_profile=squeeze(sum(sum(v)));
            
            
            
            %%  chose vertical coord and orientation for interpolation
            % according to vertical coordinates of input data and schism
            % a) Vertical coordinates are Z -> transform to sigma
            
            %             if vrt_in_sigma==0  % uniform depth levels
            %                 [sigma_selfe]=calculate_z_coordinate_selfe(bd_salt(i,4),nrvt);
            %                 zselfe=bd_salt(i,end)*sigma_selfe;
            %             else   %b) % Vertical coordinate is already sigma
            %                 [zselfe]=calculate_z_coordinate_selfe(bd_salt(i,4),nrvt);
            %             end
            
            zselfe=zsselfe(:,i);
    
            %perform vertical interpolation
            if flexible_vertical_coordinate==0 %OLD VERTICAL COORDINATE
                try
                ivalid=~isnan(salt_profile);
                salt_profile_intp=interp1(z_salt(ivalid),salt_profile(ivalid),zselfe,vert_interp_meth,'extrap');
                ivalid=~isnan(temp_profile);
                temp_profile_intp=interp1(z_temp(ivalid),temp_profile(ivalid),zselfe,vert_interp_meth,'extrap');
                ivalid=~isnan(u_profile);
                %JP quick fix: z_vel was not a profile, changed as z_temp
                u_profile_intp=interp1(z_temp(ivalid),u_profile(ivalid),zselfe,vert_interp_meth,'extrap');
                v_profile_intp=interp1(z_temp(ivalid),v_profile(ivalid),zselfe,vert_interp_meth,'extrap');
                catch err
                    %err.message
                    %keyboard
                    % use previouse profile if error
                    %elev_interp=valsarray_salt_nc(1,1:nrvt,i-1,t_total);
                    %temp_profile_intp=valsarray_temp_nc(1,1:nrvt,i-1,t_total);
                    %u_profile_intp=valsarray_vel_nc(1,1:nrvt,i-1,t_total);
                    %v_profile_intp=valsarray_vel_nc(2,1:nrvt,i-1,t_total);
                    
                end
            else %NEW VERTICAL flexibile COORDINATE
                ivalid=~isnan(salt_profile);
                salt_profile_intp(ibottom(i):nrvt)=interp1(z_salt(ivalid),salt_profile(ivalid),zselfe(ibottom(i):end),vert_interp_meth,'extrap');
                ivalid=~isnan(temp_profile);
                temp_profile_intp(ibottom(i):nrvt)=interp1(z_temp(ivalid),temp_profile(ivalid),zselfe(ibottom(i):end),vert_interp_meth,'extrap');
                ivalid=~isnan(u_profile);
                u_profile_intp(ibottom(i):nrvt)=interp1(z_temp(ivalid),u_profile(ivalid),zselfe(ibottom(i):end),vert_interp_meth,'extrap');
                v_profile_intp(ibottom(i):nrvt)=interp1(z_temp(ivalid),v_profile(ivalid),zselfe(ibottom(i):end),vert_interp_meth,'extrap');
                
                % set values in sigma layers below below bottom to bottom
                salt_profile_intp(1:ibottom(i)-1)=salt_profile_intp(ibottom(i));
                temp_profile_intp(1:ibottom(i)-1)=temp_profile_intp(ibottom(i));
                u_profile_intp(1:ibottom(i)-1)=u_profile_intp(ibottom(i));
                v_profile_intp(1:ibottom(i)-1)=v_profile_intp(ibottom(i));
            end
            
    
            
            %% assing to arrays
            %try
            valsarray_elev_nc(1,1,i,t_total)=elev_interp;
            valsarray_salt_nc(1,1:nrvt,i,t_total)=salt_profile_intp;
            valsarray_temp_nc(1,1:nrvt,i,t_total)=temp_profile_intp;
            valsarray_vel_nc(1,1:nrvt,i,t_total)=u_profile_intp;
            valsarray_vel_nc(2,1:nrvt,i,t_total)=v_profile_intp;
            %catch err
            %   err.message
            %   keyboard
            %end
            
        end  % END BOUNDARY NODE LOOP
        
    end % time
    
    netcdf.close(ncids.elev);
    netcdf.close(ncids.salt);
    netcdf.close(ncids.temp);
    netcdf.close(ncids.vel);
end %day loop


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%   Reoder and Write Variable arrays %%%%%%%%%%%%%%%%%%%%


% flup ud arrays to have surface values at higher z indices
valsarray_salt_nc=valsarray_salt_nc(1,end:-1:1,:,:);
valsarray_temp_nc=valsarray_temp_nc(1,end:-1:1,:,:);
valsarray_vel_nc=valsarray_vel_nc(:,end:-1:1,:,:);


%% interpolation to desired time step
tforce=(0:t_total-1)'*forcing_dt+offset_time;
tinterp=(0:t_total-2)'*forcing_dt;   %(tforce(1):selfe_dt:tforce(end))'; fix with zeor dt starting




% padd fpr restart
valsarray_elev_nc2=zeros(1,1,nbd,n_total);
valsarray_salt_nc2=zeros(1,nrvt,nbd,n_total);
valsarray_temp_nc2=zeros(1,nrvt,nbd,n_total);
valsarray_vel_nc2=zeros(2,nrvt,nbd,n_total);


%% interpolate linear in time
ti=1;
dt=[tinterp(ti)-tforce(ti),tforce(ti+1)-tinterp(ti)];
w=(1./dt)/sum((1./dt));
w(isnan(w))=1; if sum(w)>1;keyboard;end
for ti=1:length(tinterp)-(tforce(end)<tinterp(end));
    valsarray_elev_nc2(:,:,:,ti)=valsarray_elev_nc(:,:,:,ti)*w(1) +valsarray_elev_nc(:,:,:,ti+1)*w(2);
    valsarray_salt_nc2(:,:,:,ti)=valsarray_salt_nc(:,:,:,ti)*w(1) +valsarray_salt_nc(:,:,:,ti+1)*w(2);
    valsarray_temp_nc2(:,:,:,ti)=valsarray_temp_nc(:,:,:,ti)*w(1) +valsarray_temp_nc(:,:,:,ti+1)*w(2);
    valsarray_vel_nc2(:,:,:,ti)=valsarray_vel_nc(:,:,:,ti)*w(1) +valsarray_vel_nc(:,:,:,ti+1)*w(2);
end

if tforce(end)<tinterp(end)
    disp('delete last time step. not covered by forcing data i.e. not endint at 24:00')
    valsarray_elev_nc2(:,:,:,end)=[];
    valsarray_salt_nc2(:,:,:,end)=[];
    valsarray_temp_nc2(:,:,:,end)=[];
    valsarray_vel_nc2(:,:,:,end)=[];
    tinterp=tinterp(1:end-1);
end

%% zero padding if restart != start
ntotal=length(tinterp);
if startdate==restartdate
    time_nc=(0:ntotal-1)*selfe_dt;
else
    
    timestepadd=donedays*86400; %/selfe_dt;
    npadd=timestepadd/selfe_dt;
    time_nc=(0:ntotal-2 +  npadd )*selfe_dt +  timestepadd; %assume previous rund start at 00:00 and ends at 24:00
    size(valsarray_elev_nc2)                                                                       %threfore remove first of current time steps in arrays
    valsarray_elev_nc2(:,:,:,1)=[];
    valsarray_salt_nc2(:,:,:,1)=[];
    valsarray_temp_nc2(:,:,:,1)=[];
    valsarray_vel_nc2(:,:,:,1)=[];
    
    % zero padding at beginning
    valsarray_elev_nc2=cat(4,zeros(1,1,nbd,npadd ) ,valsarray_elev_nc2);
    valsarray_salt_nc2=cat(4,zeros(1,nrvt,nbd,npadd ) ,valsarray_salt_nc2);
    valsarray_temp_nc2=cat(4,zeros(1,nrvt,nbd,npadd ) ,valsarray_temp_nc2);
    valsarray_vel_nc2=cat(4,zeros(2,nrvt,nbd,npadd ) ,valsarray_vel_nc2);
end
ntotal=length(time_nc);


%% writing


% initialize
%output - define netcdf output %VH 16.10.2018 after JP,26/03/2018
ncids.Tout=netcdf.create(fullfile(bnddir,['TEM_3D' suffix daterange  '.th.nc']),'NC_WRITE'); %biogc tracers
%define dimensions
obn_dimID = netcdf.defDim(ncids.Tout,'nOpenBndNodes',sum(Len_force_bd));
nrvt_dimID = netcdf.defDim(ncids.Tout,'nLevels',nrvt);
comp_dimID = netcdf.defDim(ncids.Tout,'nComponents',1);
one_dimID = netcdf.defDim(ncids.Tout,'one',1);
time_dimID = netcdf.defDim(ncids.Tout,'time',netcdf.getConstant('NC_UNLIMITED'));
%define variables
Ttstep_varid = netcdf.defVar(ncids.Tout,'time_step','NC_DOUBLE',one_dimID);
netcdf.putAtt(ncids.Tout,Ttstep_varid,'long_name','time step in seconds');
Ttime_varid = netcdf.defVar(ncids.Tout,'time','NC_DOUBLE',time_dimID);
netcdf.putAtt(ncids.Tout,Ttime_varid,'long_name','simulation time in seconds');
tem_varid = netcdf.defVar(ncids.Tout,'time_series','NC_DOUBLE',[comp_dimID nrvt_dimID obn_dimID time_dimID]);
netcdf.endDef(ncids.Tout)

ncids.Sout=netcdf.create(fullfile(bnddir,['SAL_3D' suffix daterange  '.th.nc']),'NC_WRITE'); %biogc tracers
%define dimensions
obn_dimID = netcdf.defDim(ncids.Sout,'nOpenBndNodes',sum(Len_force_bd));
nrvt_dimID = netcdf.defDim(ncids.Sout,'nLevels',nrvt);
comp_dimID = netcdf.defDim(ncids.Sout,'nComponents',1);
one_dimID = netcdf.defDim(ncids.Sout,'one',1);
time_dimID = netcdf.defDim(ncids.Sout,'time',netcdf.getConstant('NC_UNLIMITED'));
%define variables
Ststep_varid = netcdf.defVar(ncids.Sout,'time_step','NC_DOUBLE',one_dimID);
netcdf.putAtt(ncids.Sout,Ststep_varid,'long_name','time step in seconds');
Stime_varid = netcdf.defVar(ncids.Sout,'time','NC_DOUBLE',time_dimID);
netcdf.putAtt(ncids.Sout,Stime_varid,'long_name','simulation time in seconds');
sal_varid = netcdf.defVar(ncids.Sout,'time_series','NC_DOUBLE',[comp_dimID nrvt_dimID obn_dimID time_dimID]);
netcdf.endDef(ncids.Sout)

ncids.Uout=netcdf.create(fullfile(bnddir,['uv3D' suffix daterange  '.th.nc']),'NC_WRITE'); %biogc tracers
%define dimensions
obn_dimID = netcdf.defDim(ncids.Uout,'nOpenBndNodes',sum(Len_force_bd));
nrvt_dimID = netcdf.defDim(ncids.Uout,'nLevels',nrvt);
comp_dimID = netcdf.defDim(ncids.Uout,'nComponents',2);
one_dimID = netcdf.defDim(ncids.Uout,'one',1);
time_dimID = netcdf.defDim(ncids.Uout,'time',netcdf.getConstant('NC_UNLIMITED'));
%define variables
Utstep_varid = netcdf.defVar(ncids.Uout,'time_step','NC_DOUBLE',one_dimID);
netcdf.putAtt(ncids.Uout,Utstep_varid,'long_name','time step in seconds');
Utime_varid = netcdf.defVar(ncids.Uout,'time','NC_DOUBLE',time_dimID);
netcdf.putAtt(ncids.Uout,Utime_varid,'long_name','simulation time in seconds');
vel_varid = netcdf.defVar(ncids.Uout,'time_series','NC_DOUBLE',[comp_dimID nrvt_dimID obn_dimID time_dimID]);
netcdf.endDef(ncids.Uout)

ncids.Eout=netcdf.create(fullfile(bnddir,['elev2D' suffix daterange  '.th.nc']),'NC_WRITE'); %biogc tracers
%define dimensions
obn_dimID = netcdf.defDim(ncids.Eout,'nOpenBndNodes',sum(Len_force_bd));
nrvt_dimID = netcdf.defDim(ncids.Eout,'nLevels',1);
comp_dimID = netcdf.defDim(ncids.Eout,'nComponents',1);
one_dimID = netcdf.defDim(ncids.Eout,'one',1);
time_dimID = netcdf.defDim(ncids.Eout,'time',netcdf.getConstant('NC_UNLIMITED'));
%define variables
Etstep_varid = netcdf.defVar(ncids.Eout,'time_step','NC_DOUBLE',one_dimID);
netcdf.putAtt(ncids.Eout,Etstep_varid,'long_name','time step in seconds');
Etime_varid = netcdf.defVar(ncids.Eout,'time','NC_DOUBLE',time_dimID);
netcdf.putAtt(ncids.Eout,Etime_varid,'long_name','simulation time in seconds');
elev_varid = netcdf.defVar(ncids.Eout,'time_series','NC_DOUBLE',[comp_dimID nrvt_dimID obn_dimID time_dimID]);
netcdf.endDef(ncids.Eout)




% write array output
netcdf.putVar(ncids.Eout,Etstep_varid, selfe_dt)         % VH
netcdf.putVar(ncids.Eout,Etime_varid, 0,ntotal, time_nc)
netcdf.putVar(ncids.Eout,elev_varid, [0 0 0 0],[1 1 sum(Len_force_bd) ntotal+1], valsarray_elev_nc2)

netcdf.putVar(ncids.Uout,Utstep_varid, selfe_dt)
netcdf.putVar(ncids.Uout,Utime_varid, 0,ntotal, time_nc)
netcdf.putVar(ncids.Uout,vel_varid, [0 0 0 0],[2 nrvt sum(Len_force_bd) ntotal+1],  valsarray_vel_nc2)

netcdf.putVar(ncids.Sout,Stime_varid, 0,ntotal, time_nc)
netcdf.putVar(ncids.Sout,Ststep_varid, selfe_dt)
netcdf.putVar(ncids.Sout,sal_varid, [0 0 0 0],[1 nrvt sum(Len_force_bd) ntotal+1], valsarray_salt_nc2)

netcdf.putVar(ncids.Tout,Ttstep_varid, selfe_dt)
netcdf.putVar(ncids.Tout,Ttime_varid, 0,ntotal, time_nc)
netcdf.putVar(ncids.Tout,tem_varid, [0 0 0 0],[1 nrvt sum(Len_force_bd) ntotal+1],  valsarray_temp_nc2)


netcdf.close(ncids.Eout)
netcdf.close(ncids.Sout)
netcdf.close(ncids.Tout)
netcdf.close(ncids.Uout)


fclose all;
cd(bnddir)
toc


end
%% Hilfsfunktionen

function varid = get_varid_from_standard_name(ncid,standardname)
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
varid=-1;
for i = 0:nvars-1
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncid,i);
    if strcmp(netcdf.getAtt(ncid,i,'standard_name'),standardname);
        varid=i;
        break
    end
end
if varid==-1;disp('varname not found');keyboard; error('varname not found');end
end

function [sigma_selfe]=calculate_z_coordinate_selfe(depth_selfe,number_levels)


sigma_level=linspace(0,-1,number_levels);

h_c=10;

theta_b=0;

theta_f=4;

for sig=1:number_levels
    
    cs(sig)=(1-theta_b)*sinh(theta_f*sigma_level(sig))/sinh(theta_f)+theta_b*(tanh(theta_f*(sigma_level(sig)+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5);
    
end

if depth_selfe<=h_c
    
    z_selfe=sigma_level.*depth_selfe;
    
else
    
    z_selfe=h_c.*sigma_level+(depth_selfe-h_c).*cs;
    
end

sigma_selfe=0.+cumsum(diff(z_selfe)./depth_selfe);
sigma_selfe=[0 sigma_selfe];

return

end





% Translate profiles into sigma frame
function [ sigma_bnd_profile sigma_selfe bnd_node_profiles_temp]=profile_Z2sigma(bd,bnd_node_profiles,elev_interp,i,depth,nz,nrvt)

%keyboard
bnd_node_profiles_temp=bnd_node_profiles(:,i);
depth_temp=depth;


ind_nan_vert=isnan(bnd_node_profiles(:,i));

if sum(ind_nan_vert)==nz %alles nans
    % %         i
    % %         bnd_node_profiles_profileIn_new(:,i)=bnd_node_profiles_profileIn_new(:,i-1);
    % %         switch_nan=1;
    % %         bd_profileIn(i,4)
    % %         bd_profileIn(i-1,4)
    
    ind_nans=1;
    
    %nehme profil von vorhergehenden punkt solange profil mit
    %nicht nan werten erreicht
    while sum(ind_nan_vert)==nz
        clear bnd_node_profiles_temp
        clear depth_temp
        clear ind_nan_vert
        bnd_node_profiles_temp=bnd_node_profiles(:,i-ind_nans);
        depth_temp=depth;
        ind_nan_vert=isnan(bnd_node_profiles_temp);
        %i
        ind_nans=ind_nans+1;
        
    end
    
end

%Reduziere Wassers???ule auf maximale Tiefe der benachbarten
%s???ule
bnd_node_profiles_temp(ind_nan_vert)=[];
depth_temp(ind_nan_vert)=[];
depth_old=depth_temp;
bnd_node_profiles_old=bnd_node_profiles_temp;
[~,ind_dist]=min(abs(bd(i,4)-depth_temp));
if depth_temp(ind_dist)>bd(i,4)
    ind_dep=depth_temp>=depth_temp(ind_dist);
    depth_temp(ind_dep)=[];
    bnd_node_profiles_temp(ind_dep)=[];
    depth_temp(ind_dist)=bd(i,4);
    bnd_node_profiles_temp(ind_dist)=interp1(depth_old,bnd_node_profiles_old,depth_temp(ind_dist),'nearest','extrap');   %bnd_node_profiles_temp(ind_dist-1);
elseif depth_temp(ind_dist)<bd(i,4)   %Boundary at NN deeper than Forcing grid
    ind_dep=depth_temp>depth_temp(ind_dist);
    depth_temp(ind_dep)=[];
    bnd_node_profiles_temp(ind_dep)=[];
    depth_temp(ind_dist+1)=bd(i,4);
    bnd_node_profiles_temp(ind_dist+1)=interp1(depth_old,bnd_node_profiles_old,depth_temp(ind_dist+1),'nearest','extrap'); %bnd_node_profiles_temp(ind_dist);
elseif depth_temp(ind_dist)==bd(i,4)
    ind_dep=depth_temp>depth_temp(ind_dist);
    depth_temp(ind_dep)=[];
    bnd_node_profiles_temp(ind_dep)=[];
    %depth_temp(ind_dist)=bd(i,4);
    %bnd_node_profiles_temp(ind_dist)=interp1(depth_old,bnd_node_profiles_old,depth_temp(ind_dist),'nearest','extrap');   %bnd_node_profiles_temp(ind_dist-1);
end

sigma_bnd_profile=-1*(depth_temp+elev_interp(i))./(max(depth_temp)+elev_interp(i));

[sigma_selfe]=calculate_z_coordinate_selfe(bd(i,4),nrvt);



end


function [x,y]=ll2xy(lat,lon,reflat,reflong)  % from batry
%ll2xy::convert latitude and longitude to distance in meters
%
% [X,Y]=LL2XYGOM(LAT,LON,REFLAT,REFLON)
%
% Returns X and Y vectors (as distance from arbitrary
% origin REFLAT and REFLON) for vectors of latitudes
% and longitudes, LAT and LON.
%
% REFLAT and REFLON default to Boston
%
% LAT and LON may be specified as LON+i*LAT
%
% Specifying only a single output yields X+i*Y

% CVL, 7-10-97
% Hacked from mercgom2 from C. Naimie.

r=6.3675E+6;

if nargin<3
    reflong=-71.03*pi/180;
    reflat=42.35*pi/180;
end

if nargin==4,%convert to radians
    reflong=reflong*pi/180;
    reflat=reflat*pi/180;
end

if nargin==1
    lon=real(lat);
    lat=imag(lat);
end

xo=r*cos(reflat)*reflong;
yo=r*cos(reflat)*log((1.0+sin(reflat))/cos(reflat));

rlong=lon*pi/180;
rlat=lat*pi/180;
x=r*cos(reflat).*rlong-xo;
y=r*cos(reflat).*log((1.0+sin(rlat))./cos(rlat))-yo;

if nargout<=1
    x=x+i*y;
end
end
