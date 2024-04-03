"""
Quick-and-dirty HVC
"""
from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.units as u
from astropy.io.votable import parse_single_table
import subprocess
import matplotlib.pyplot as plt
import numpy as np

"""
This script simply consists of steps to cobbled together. It is not to just be run start to finish
"""

"""
Filter Moss+2013 table to objects of a certain area
"""

def format_sexagesimal(coord):
    ra, dec = coord.split()[:3], coord.split()[3:]
    ra_formatted = ':'.join(ra)
    dec_formatted = f"{dec[0]}:{dec[1]}:{dec[2]}"
    return f"{ra_formatted} {dec_formatted}"

# Read the VOTable file
table = parse_single_table('/Users/canders/Projects/HVCs/catalogues/vizier_Moss2013_HVCs.vot').to_table()

# Mask to decent sized area
mask = (table['Area']>1) & (table['Area']<np.pi)

table_big_HVCs = table[mask]

# Assuming the coordinates are in columns 'ra' and 'dec', adjust as needed
#ra_dec = [f"{row['RAJ2000']} {row['DEJ2000']}" for row in table_big_HVCs]
ra_dec = [format_sexagesimal(f"{row['RAJ2000']} {row['DEJ2000']}") for row in table_big_HVCs]

# Write out
with open('/Users/canders/Projects/HVCs/big_HVC_coordinates_area_1-1p5_sqdeg.txt', 'w') as file:
    for coordinate in ra_dec:
        file.write(coordinate + '\n')

"""
Now run the POSSUM X-matcher
"""
result = subprocess.run(['python', '/Users/canders/Dropbox/Scripts/python/POSSUM_coverage_cmdline.py', '/Users/canders/Projects/HVCs/big_HVC_coordinates_area_1-1p5_sqdeg.txt'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
return_code = result.returncode
stdout = result.stdout
stderr = result.stderr

"""
Now load in the output from the X-matcher 
"""
observed_HVCs_table = Table.read('/Users/canders/Projects/HVCs/Observed_HVC_positions_SBIDs_beams.txt',format='ascii')


"""
Here, you need to run in a terminal where the table "t" has been defined from my POSSUM groups project. Need to amend this by implementing a load-in feature.
"""


"""
Calculate RM values and offsets from HVCs.
"""

# Arrays to store the results
l_offsets = []
m_offsets = []
separations = []
faraday_depths = []

# Loop through HVCs
for hvc in observed_HVCs_table:
    
    print('Examining object...')
    
    ra_deg = hvc['ra_deg']
    dec_deg = hvc['dec_deg']
    nearest_SBID = hvc['nearest_SBID']
    
    if str(nearest_SBID) not in list(set(t['sbid'])):
    	print('	Not in t (t needs updating with latest obs)...')
    	continue
    
    # Mask "t" based on the nearest_SBID value of the HVC
    mask = (t['sbid']==str(nearest_SBID)) & (t['LMC_sep_deg']>30)
    
    filtered_t = t[mask]
    
    # Define HVC coordinate
    hvc_coord = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg)
    
    # Loop through filtered radio RM_sources
    for RM_source in filtered_t:
        RM_source_coord = RM_source['ra_dec_obj']
        
        # Calculate on-sky separation
        separation = hvc_coord.separation(RM_source_coord).to('deg')
        if separation < 6*u.deg:
            # Calculate l and m direction offsets
            l, m = hvc_coord.spherical_offsets_to(RM_source_coord)
            l_offsets.append(l.to('deg').value)
            m_offsets.append(m.to('deg').value)
            separations.append(separation.value)
            faraday_depths.append(RM_source['faraday_depth_radmm'])
            
# Create a scatter plot with size proportional to 'faraday_depth_radmm'
colors = ['red' if depth > 0 else 'blue' for depth in faraday_depths]
plt.scatter(l_offsets, m_offsets, c=colors, s=[abs(depth) for depth in faraday_depths])
plt.xlabel('l offset (degrees)')
plt.ylabel('m offset (degrees)')
pl.axhline(y=0,c='r')
pl.axvline(x=0,c='r')
plt.show()