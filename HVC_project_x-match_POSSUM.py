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
#mask = (table['Area']>1) & (table['Area']<np.pi)
mask = table['Area']>1

table_big_HVCs = table[mask]

# Assuming the coordinates are in columns 'ra' and 'dec', adjust as needed
#ra_dec = [f"{row['RAJ2000']} {row['DEJ2000']}" for row in table_big_HVCs]
ra_dec = [format_sexagesimal(f"{row['RAJ2000']} {row['DEJ2000']}") for row in table_big_HVCs]

# Write out
with open('/Users/canders/Projects/HVCs/big_HVC_coordinates_area_greater_than_1_sqdeg.txt', 'w') as file:
    for coordinate in ra_dec:
        file.write(coordinate + '\n')

"""
Now run the POSSUM X-matcher
"""
result = subprocess.run(['python', '/Users/canders/Dropbox/Scripts/python/POSSUM_coverage_cmdline.py', '/Users/canders/Projects/HVCs/big_HVC_coordinates_area_greater_than_1_sqdeg.txt'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
return_code = result.returncode
stdout = result.stdout
stderr = result.stderr

"""
Now load in the output from the X-matcher 
"""
observed_HVCs_table = Table.read('/Users/canders/Projects/HVCs/Observed_HVC_positions_SBIDs_beams.txt',format='ascii')