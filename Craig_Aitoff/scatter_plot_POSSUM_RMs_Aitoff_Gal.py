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

# Data from astropy table
rm_centres = np.array(combined_table['ra_dec_deg'])
ra_deg = rm_centres[:, 0]
dec_deg = rm_centres[:, 1]
faraday_depth_radmm = combined_table['RM']

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
norm = Normalize(vmin=-100, vmax=100)

# Plot in Aitoff projection
fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection': 'mollweide'})
sc = ax.scatter(l, b, c=faraday_depth_radmm, marker='.', cmap='plasma', norm=norm, s=2, alpha=0.6, edgecolors='none')

# Colorbar at the bottom, horizontal and thinner
sm = plt.cm.ScalarMappable(cmap='RdBu_r', norm=norm)
#sm = plt.cm.ScalarMappable(cmap='plasma', norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, orientation='horizontal', pad=0.05, aspect=50)
cbar.set_label('Faraday Depth (rad/m^2)', fontsize=12)

# Labels and title
ax.set_xlabel('Galactic Longitude (l)')
ax.set_ylabel('Galactic Latitude (b)')
ax.set_title('POSSUM RMs @ 26 May 2024')

# Make background light gray
ax.set_facecolor('lightgray')

# Convert to hours
ax.xaxis.set_major_formatter(plt.FuncFormatter(galactic_to_hours))

# Activate faint grid lines
ax.grid(True, which='both', color='gray', linestyle='--', linewidth=0.5, alpha=0.5)

# Save the figure as a PDF with close borders
plt.savefig('POSSUM_RMs_Aitoff_26May_Moll_gal_FG_corr_abs_zorder.pdf', format='pdf', bbox_inches='tight')

plt.show()
