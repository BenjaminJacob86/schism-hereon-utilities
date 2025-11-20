
import sys
from schism import *
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
plt.ion()
s=schism_setup()
#ph,ch,ax=s.plotAtnodesGeo(np.asarray(s.depths),region_limit=[28.7,29.4,40.9,41.35],extend=None,land_color='white')

plt.rcParams.update({'font.size': 11})


#plt.tripcolor(s.x,s.y,s.nvplt,facecolors=s.depths,shading='gouraud',cmap=plt.cm.jet)# shading needs gouraud to allow correct update

#gl.ylabel_style = {"size": 11}

ph,ch,ax=s.plotAtnodesGeo(np.asarray(s.depths),proj=ccrs.Mercator(),region_limit=[28.7,29.4,40.9,41.35],extend=None,landcolor='white',cmap=plt.cm.jet)

#,proj=ccrs.PlateCarree()

ch.set_label('Depth [m]')
ph.set_clim((0,100))

ch.set_ticks(np.arange(0,110,10))
ch.set_ticklabels(np.arange(0,110,10))


nxticks,nyticks = 8,6

# Region in lon/lat (same as your region_limit)
lon_min, lon_max, lat_min, lat_max = [28.7, 29.4, 40.9, 41.35]

# Generate evenly spaced longitude and latitude ticks
xticks = np.linspace(lon_min, lon_max, nxticks)
yticks = np.linspace(lat_min, lat_max, nyticks)

# Tell Cartopy to place ticks in PlateCarree (lon/lat),
# while transforming them into your UTM projection
ax.set_xticks(xticks, crs=ccrs.PlateCarree())
ax.set_yticks(yticks, crs=ccrs.PlateCarree())

# Format as longitude/latitude
ax.xaxis.set_major_formatter(LongitudeFormatter(number_format=".1f"))
ax.yaxis.set_major_formatter(LatitudeFormatter(number_format=".1f"))


plt.tight_layout()
plt.savefig('Bosporus2.png',dpi=300)

# Optional styling
ax.tick_params(labelsize=10)













nxticks,nyticks = 8,8

zoom_extend=[28.7,29.4,40.9,41.35]
xticks=np.unique(np.round(np.linspace((zoom_extend[0]),(zoom_extend[1]),nxticks),1))
yticks=np.unique(np.round(np.linspace((zoom_extend[2]),(zoom_extend[3]),nyticks),1))


proj=ax.projection
ax.set_xticks(xticks, crs=proj)
ax.set_yticks(yticks, crs=proj)
lon_formatter = LongitudeFormatter(number_format='.1f',degree_symbol='',dateline_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.1f',degree_symbol='')
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)