function combine_sed_TS

% first, create a hotstart file contans T, S with AMM15 data (use genHot_AMM15.py )
% second, create a hotstart file with sediment (use ../script/createLastHot.batch )
% third, run this function to combine AMM15 TS with sediment
% (specify input_dir, input_ts, input_sed)


input_dir = '/work/gg0028/g260114/RUNS/GermanBight/GB_timeslices/GB1982_wave/split/GB_dt80';
input_ts  = 'hotstart.nc_hydro';
input_sed  = 'hotstart.nc_sed_cold_dt1_test';
input_sed  = 'hotstart.nc_sed_cold_dt1';
tr_nd_ts   = squeeze(double(ncread([input_dir, '/', input_ts],'tr_nd')));
tr_nd0_ts  = squeeze(double(ncread([input_dir, '/', input_ts],'tr_nd0')));
tr_el_ts   = squeeze(double(ncread([input_dir, '/', input_ts],'tr_el')));

%tr_nd   = squeeze(double(ncread([input_dir, '/hotstarts/', input_sed],'tr_nd')));
%tr_nd0   = squeeze(double(ncread([input_dir, '/hotstarts/', input_sed],'tr_nd0')));
%tr_el   = squeeze(double(ncread([input_dir, '/hotstarts/', input_sed],'tr_el')));

tr_nd   = squeeze(double(ncread([input_dir,  '/',input_sed],'tr_nd')));
tr_nd0   = squeeze(double(ncread([input_dir, '/', input_sed],'tr_nd0')));
tr_el   = squeeze(double(ncread([input_dir,  '/',input_sed],'tr_el')));



tr_nd(1:2,:,:) = tr_nd_ts;
tr_nd0(1:2,:,:)= tr_nd0_ts;
tr_el(1:2,:,:)= tr_el_ts;


%ncwrite([input_dir,'/hotstarts/',input_sed],'tr_nd',tr_nd);
%ncwrite([input_dir,'/hotstarts/',input_sed],'tr_nd0',tr_nd0);
%ncwrite([input_dir,'/hotstarts/',input_sed],'tr_el',tr_el);

ncwrite([input_dir, '/',input_sed],'tr_nd',tr_nd);
ncwrite([input_dir, '/',input_sed],'tr_nd0',tr_nd0);
ncwrite([input_dir, '/',input_sed],'tr_el',tr_el);


exit
