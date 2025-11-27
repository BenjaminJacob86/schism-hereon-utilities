import cdsapi
#import certifi
#import urllib3
import datetime as dt

# download complete month of era 5
startdate=dt.datetime(2024,6,1)
enddate=dt.datetime(2024,6,4)


variables=['10m_u_component_of_wind', '10m_v_component_of_wind', '2m_dewpoint_temperature', '2m_temperature', 'mean_sea_level_pressure', 'sea_surface_temperature', 'snowfall', 'surface_net_thermal_radiation', 'surface_solar_radiation_downwards', 'surface_thermal_radiation_downwards', 'total_cloud_cover', 'total_precipitation', 'sea_ice_cover']

#variables=['10m_u_component_of_wind', '10m_v_component_of_wind']
product='reanalysis-era5-single-levels'

days=['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31']
#days=['01']
#product='reanalysis-era5-pressure-levels'
#variables=['specific_humidity']


year=startdate.year
mon=startdate.month
date=dt.datetime(year,mon,1)
while date <= enddate:
	c = cdsapi.Client()
	print('downloading ' + product+ '' +str(date))
	if' pressure' in product:
		c.retrieve(
		product,
		{
			'product_type':'reanalysis',
			'format':'netcdf',
			'variable':variables,
			'pressure_level':'1000',
			'grid': "0.25/0.25",
	        	'area': '81/-120/-10/68',
	        	'year':'{:d}'.format(year),
		        'month':'{:02d}'.format(mon),
        		'day':days,
		        'time':['00:00', '01:00', '02:00', '03:00', '04:00', '05:00', '06:00', '07:00', '08:00', '09:00', '10:00', '11:00', '12:00', '13:00', '14:00', '15:00', '16:00', '17:00', '18:00', '19:00', '20:00', '21:00', '22:00', '23:00']
        	},
	        'ERA5CDS_spfh{:4d}{:02d}.nc'.format(year,mon))
	else:
		c.retrieve(
		product,
		{
			'product_type':'reanalysis',
			'format':'netcdf',
			'variable':variables,
			'grid': "0.25/0.25",
	        	'area': '81/-120/-10/68',
	        	'year':'{:d}'.format(year),
		        'month':'{:02d}'.format(mon),
        		'day':days,
		        'time':['00:00', '01:00', '02:00', '03:00', '04:00', '05:00', '06:00', '07:00', '08:00', '09:00', '10:00', '11:00', '12:00', '13:00', '14:00', '15:00', '16:00', '17:00', '18:00', '19:00', '20:00', '21:00', '22:00', '23:00']
        	},
	        'ERA5CDS{:4d}{:02d}.nc'.format(year,mon))


	if mon < 12:
		mon+=1
	else:
		year+=1
		mon=1
	date=dt.datetime(year,mon,1)

