import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import t as student_t
#import seaborn as sns

# === SETUP === #

# Get RM column

from astropy.io import ascii

def get_RMs(file, param):
    data = ascii.read(file)
    col = data["RM"].data - data[param].data
    col = col[~np.isnan(col)]
    return col

rm_column = get_RMs("../data_processed/proc_rms.ecsv", "interpolation_raw")
rm_column2 = get_RMs("../data_processed/proc_rms.ecsv", "interpolation_cor")

# Set style using seaborn
#sns.set_style("whitegrid")
#palette = sns.color_palette("icefire")

# === MATPLOTLIB SETUP === #

# Set figure size suitable for single column in a two-column paper
plt.figure(figsize=(8.5, 7))
ax = plt.gca()

# Make spines less prominent
#ax.spines['top'].set_linewidth(0.3)
#ax.spines['right'].set_linewidth(0.3)
#ax.spines['bottom'].set_linewidth(0.3)
#ax.spines['left'].set_linewidth(0.3)

# Add minor ticks
#ax.minorticks_on()

# Set the background color to white
ax.set_facecolor('white')

# Set grid properties and ensure it appears behind all other elements
#ax.grid(color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
#ax.set_axisbelow(True)

# === BOOTSTRAP & T-DISTRIBUTION === #

num_samples = 10

def bootstrap_fit(data, num_bootstrap_samples):
    params_bootstrap = np.zeros((num_bootstrap_samples, 3))
    
    for i in range(num_bootstrap_samples):
        print(f'Sample {i} of {num_bootstrap_samples}')
        bootstrap_sample = np.random.choice(data, len(data), replace=True)
        params_bootstrap[i, :] = student_t.fit(bootstrap_sample)
    
    params_mean = params_bootstrap.mean(axis=0)
    params_std = params_bootstrap.std(axis=0)
    
    return params_mean, params_std

# Fit the student's t-distribution to your data and plot the PDF
def bootstrap_to_t(column, name, show=True):
    print("Fitting t-distribution for: "+name)
    params_mean, params_std = bootstrap_fit(column, num_samples)
    df_mean, loc_mean, scale_mean = params_mean
    df_std, loc_std, scale_std = params_std

    if show:
        x = np.linspace(-100, 100, 1000)
        pdf = student_t.pdf(x, df_mean, loc_mean, scale_mean)
        ax.plot(x, pdf, label=name, linewidth=1)

    return params_mean, params_std

params_mean1, params_std1 = bootstrap_to_t(rm_column, "Raw Interpolation Residuals t-dist")
params_mean2, params_std2 = bootstrap_to_t(rm_column2, "Crosshatch-Bandpassed Residuals t-dist")

# === HISTOGRAM CREATION === #

# Create a histogram
bins = np.linspace(-100, 100, 100)
n, bins, patches = plt.hist(rm_column, bins=bins, density=True, alpha=0.5, edgecolor='none', color='#7cbee1', label="Raw Interpolation Residuals")
n2, bins2, patches2 = plt.hist(rm_column2, bins=bins, density=True, alpha=0.5, edgecolor='none', color='#e19f7c', label="Crosshatch-Bandpassed Residuals")

# Adding x and y axis labels
plt.xlabel(r'RM Residuals [rad m$^{-2}$]')
plt.ylabel('Density')

# Set x-axis limits
plt.xlim(-75, 75)

# === T-DISTRIBUTION OUTPUT === #

# Add fit parameters in a text box above the axes
def make_txt(params_mean, params_std, label):
    df_mean, loc_mean, scale_mean = params_mean
    df_std, loc_std, scale_std = params_std
    #textstr = '\n'.join((
    #    r'$\mathrm{df}=%.2f \pm %.2f$' % (df_mean, df_std),
    #    r'$\mathrm{loc}=%.2f \pm %.2f \, \mathrm{[rad \, m^{-2}]}$' % (loc_mean, loc_std),
    #    r'$\mathrm{scale}=%.2f \pm %.2f \, \mathrm{[rad \, m^{-2}]}$' % (scale_mean, scale_std)
    #))
    textstr = '\n'.join((
        str(label.upper()),
        "|--------------------------------------|",
        "| Parameter          | Mean   | Stdev  |",
        "|--------------------------------------|",
        #| Degrees of Freedom | 0.0000 | 0.0000 |
        '| Degrees of Freedom | %.4f | %.4f |' % (df_mean, df_std),
        '| Centerpoint        | %.4f | %.4f |' % (loc_mean, loc_std),
        '| Internal Scale     | %.3f | %.4f |' % (scale_mean, scale_std),
        "|--------------------------------------|"
    ))
    return textstr
    

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
textstr = make_txt(params_mean1, params_std1, "Raw Interpolation Residuals t-dist") + "\n\n" + make_txt(params_mean2, params_std2, "Crosshatch-Bandpassed Residuals t-dist")
print(textstr)
#ax.text(0.48, 1.05, textstr, transform=ax.transAxes, fontsize=7, verticalalignment='bottom', bbox=props)

# Add legend above the axes
#legend = ax.legend(["t-distribution fit", "RRM", "RM"], loc='lower left', bbox_to_anchor=(0.0, 1.05))
legend = ax.legend()

plt.tight_layout()
plt.savefig("./output.png", pad_inches=0.0)

# === CHI-SQUARED TEST WITH T-DISTRIBUTION === #

def t_dist_values(hist, bins, params):
    df, loc, scale = params
    distarr = []
    for i in range(len(hist)):
        distarr.append(student_t.pdf((bins[i+1]+bins[i])/2, df, loc=loc, scale=scale))
    return np.array(distarr)

from scipy.stats import chisquare as chisq

def cs(n, y, ddof=2):
    csq = chisq(n, np.sum(n)/np.sum(y) * y, ddof=ddof)
    return [csq[0], 1 - csq[1]]

from scipy.stats import chi2

def cs2(n, y, ddof=3):
    stat = np.sum((n-y)**2 / y)
    pv = chi2.cdf(stat, len(n)-ddof)
    return [stat, 1-pv]

thist1 = t_dist_values(n, bins, params_mean1)
chi1 = cs2(n, thist1, ddof=2)

thist2 = t_dist_values(n2, bins2, params_mean2)
chi2 = cs2(n2, thist2, ddof=2)

def make_txt_chisq(chisq, label):
    textstr = '\n'.join((
        str(label.upper()),
        "----------------------",
        #| Statistic | 0.0000 |
        '| Statistic | %.4g |' % (chisq[0]),
        '| p-value   | %.4g |' % (chisq[1]),
        "----------------------"
    ))
    return textstr

textstr = make_txt_chisq(chi1, "Raw Interpolation Residuals Chi-Squared") + "\n\n" + make_txt_chisq(chi2, "Crosshatch-Bandpassed Residuals Chi-Squared")
print()
print(textstr)