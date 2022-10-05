import numpy as np
from schism import *

# get setup object, skip reading vgrid_file
s = schism_setup(vgrid_file='')

# get list of sources (e.g. rivers)
coords = {}

coords['Elbe'] = (9.809,53.546)
coords['Weser'] = (8.497,53.390)
coords['Ems'] = (7.403,53.241)
coords['Humber'] = (-0.262,53.719)
coords['Thames'] = (0.277,51.46)

vsource = {}

vsource['Elbe'] = 716.2 
vsource['Weser'] = 327.3
vsource['Ems'] = 83.0
vsource['Humber'] = 10.0
vsource['Thames'] = 110.1

# initialize element search tree
print('  initialize element tree')
s.init_element_tree()
print('  find rivers...')

f = open('source_sink.in','w')
f.write('%d   ! number of elements with sources\n'%len(coords))

for river in coords.keys():
  lon,lat = coords[river]
  el_id = s.find_nearest_element(lon,lat)
  el_depth = s.element_depth[el_id]
  f.write('%d    ! river %s, element depth: %0.1f m\n'%(el_id,river,el_depth))

f.write('\n')
f.write('0   ! number of elements with sinks\n')
f.close()
