import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import t as student_t
import seaborn as sns

num_bootstrap_samples = 1000

def bootstrap_fit(data, num_bootstrap_samples):
    params_bootstrap = np.zeros((num_bootstrap_samples, 3))
    
    for i in range(num_bootstrap_samples):
        print(f'Sample {i} of {num_bootstrap_samples}')
        bootstrap_sample = np.random.choice(data, len(data), replace=True)
        params_bootstrap[i, :] = student_t.fit(bootstrap_sample)
    
    params_mean = params_bootstrap.mean(axis=0)
    params_std = params_bootstrap.std(axis=0)
    
    return params_mean, params_std

# Get RM column

from astropy.io import ascii

data = ascii.read("../data_processed/proc_rms.ecsv")
rm_column = data["RM"].data - data["interpolation_raw"].data #[insert RM array here]
rm_column = rm_column[~np.isnan(rm_column)]

# Set style using seaborn
sns.set_style("whitegrid")
palette = sns.color_palette("icefire")

# Set figure size suitable for single column in a two-column paper
plt.figure(figsize=(3.5, 2.8))
ax = plt.gca()

# Make spines less prominent
ax.spines['top'].set_linewidth(0.3)
ax.spines['right'].set_linewidth(0.3)
ax.spines['bottom'].set_linewidth(0.3)
ax.spines['left'].set_linewidth(0.3)

# Add minor ticks
ax.minorticks_on()

# Set the background color to white
ax.set_facecolor('white')

# Set grid properties and ensure it appears behind all other elements
ax.grid(color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
ax.set_axisbelow(True)

# Create a histogram
bins = np.linspace(-100, 100, 100)
n, bins, patches = plt.hist(rm_column, bins=bins, density=True, alpha=0.5, edgecolor='none', color='#e19f7c')
#n2, bins2, patches2 = plt.hist(rm_column2, bins=bins, density=True, alpha=0.5, edgecolor='none', color='#7cbee1')

# Fit the student's t-distribution to your data and plot the PDF
params_mean, params_std = bootstrap_fit(rm_column, 1000)
df_mean, loc_mean, scale_mean = params_mean
df_std, loc_std, scale_std = params_std
x = np.linspace(-100, 100, 1000)
pdf = student_t.pdf(x, df_mean, loc_mean, scale_mean)
ax.plot(x, pdf, label="Fit", linewidth=1, color=palette[4])

# Adding x and y axis labels
plt.xlabel('RRM [rad m$^{-2}$]')
plt.ylabel('Density')

# Set x-axis limits
plt.xlim(-75, 75)

# Add fit parameters in a text box above the axes
textstr = '\n'.join((
    r'$\mathrm{df}=%.2f \pm %.2f$' % (df_mean, df_std),
    r'$\mathrm{loc}=%.2f \pm %.2f \, \mathrm{[rad \, m^{-2}]}$' % (loc_mean, loc_std),
    r'$\mathrm{scale}=%.2f \pm %.2f \, \mathrm{[rad \, m^{-2}]}$' % (scale_mean, scale_std),
))
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.48, 1.05, textstr, transform=ax.transAxes, fontsize=7,
        verticalalignment='bottom', bbox=props)

# Add legend above the axes
legend = ax.legend(["t-distribution fit", "RRM", "RM"], loc='lower left', bbox_to_anchor=(0.0, 1.05))