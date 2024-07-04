import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scipy.spatial.distance import cdist
from astropy.coordinates import SkyCoord
import astropy.units as u

### Functions

# Convert Galactic longitude tick labels to hours
def galactic_to_hours(x, pos):
    hours = ((-x) * 12 / np.pi) % 24
    if hours in [0, 6, 12, 18]:
        return f'{int(hours)}h'
    return ''

### Params
figsize = [13,  8.4]

# Load collated data

import sys
sys.path.append('C://Users/olive/OneDrive - Australian National University/Honours-Olivia/Programs/honours/standard_modules')
sys.path.append('C://Users/olive/OneDrive - Australian National University/Honours-Olivia/Programs/honours/project')

from collation import collator

import os
os.chdir("./Craig_Aitoff")

collated_data = collator.data_whole_sky(False, load_data=["../data_processed/proc_rms_craig_annulus","../data_processed/proc_hvcs"], h1_img="../data_catalog/hi4pi-hvc-nhi-car.fits")

combined_table = collated_data["RMs"]

# Data from astropy table
rm_centres = np.array(combined_table['ra_dec_deg'])
ra_deg = rm_centres[:, 0]
dec_deg = rm_centres[:, 1]
faraday_depth_radmm = combined_table['interpolation_cor']

# Convert RA/Dec to Galactic coordinates
coords = SkyCoord(ra=ra_deg*u.degree, dec=dec_deg*u.degree, frame='icrs')
galactic_coords = coords.galactic
l = galactic_coords.l.wrap_at(180 * u.degree).radian
b = galactic_coords.b.radian

# Flip the longitude values for the x-axis to increase eastwards
l = -l

# Sort the data by faraday_depth_radmm to control z-order
sorted_indices = np.argsort(faraday_depth_radmm)
l_sorted = l[sorted_indices]
b_sorted = b[sorted_indices]
faraday_depth_radmm_sorted = faraday_depth_radmm[sorted_indices]


# Normalize Faraday depths for color mapping
norm = Normalize(vmin=-400, vmax=400)

# Plot in Aitoff projection
colormap = plt.colormaps["seismic"]
fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': 'mollweide'})
sc = ax.scatter(l, b, c=faraday_depth_radmm, marker='.', cmap=colormap, s=2, alpha=0.6, edgecolors='none', norm=norm)

# Colorbar at the bottom, horizontal and thinner
#sm = plt.cm.ScalarMappable(cmap='RdBu_r', norm=norm)
#sm = plt.cm.ScalarMappable(cmap='plasma', norm=norm)
#sm.set_array([])
cbar = plt.colorbar(sc, orientation='horizontal', pad=0.05, aspect=50)
cbar.set_label(r'Faraday Depth [$rad m^{-2}$]', fontsize=12)

# Labels and title
ax.set_xlabel('Galactic Longitude (l)')
ax.set_ylabel('Galactic Latitude (b)')
#ax.set_title('POSSUM RMs @ 26 May 2024')

# Make background light gray
ax.set_facecolor('lightgray')

# Convert to hours
ax.xaxis.set_major_formatter(plt.FuncFormatter(galactic_to_hours))

# Activate faint grid lines
ax.grid(True, which='both', color='gray', linestyle='--', linewidth=0.5, alpha=0.5)

# Save the figure as a PDF with close borders
plt.savefig('../../../Resources/Figures/POSSUM_RMs_Aitoff_Craig_Annulus.jpg', format='jpg', bbox_inches='tight')

#plt.show()
