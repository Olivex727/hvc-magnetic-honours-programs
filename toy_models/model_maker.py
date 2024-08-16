import numpy as np
import copy

# === PART 1 - CREATE VECTOR FIELD === #

def B_sim_list(xs, ys, zs, B=1):
    s = xs.shape
    xss, yss, zss = copy.deepcopy(xs), copy.deepcopy(ys), copy.deepcopy(zs)
    for xi in range(s[0]):
        for yi in range(s[1]):
            for zi in range(s[2]):
                #[zi][yi][xi]
                x = xs[xi][yi][zi]
                y = ys[xi][yi][zi]
                z = zs[xi][yi][zi]
                xv, yv, zv = B_sim(x, y, z, B)
                xss[xi][yi][zi] = xv
                yss[xi][yi][zi] = yv
                zss[xi][yi][zi] = zv
    return [xss, yss, zss]

def B_sim_list_2d(xs, y, zs, B=1):
    s = xs.shape
    xss, zss = copy.deepcopy(xs), copy.deepcopy(zs)
    for xi in range(s[0]):
        for zi in range(s[1]):
            #[zi][yi][xi]
            x = xs[xi][zi]
            z = zs[xi][zi]
            xv, yv, zv = B_sim(x, y, z, B)
            xss[xi][zi] = xv
            zss[xi][zi] = zv
    return [xss, zss]


def B_sim(x, y, z, B=1):
    r = np.sqrt(x**2+z**2)
    y0 = 1 - np.abs(y)
    zr = 1 - 0.5 * np.abs(z)
    vec = np.array([0, 0, 0])
    scale = 0
    if (np.abs(y) < 1) and (r <= 2 - 2 * np.abs(y)):
        if z < 0:
            scale = r/y0**2 * (1-np.abs(r/y0 - 1))
            return [B * scale * -z, 0, B * scale * x]
        elif z >= 0:
            scale = zr * np.abs(x)/y0**2 * (1-np.abs(np.abs(x)/y0 - 1))
            return [0, 0, B * scale * x]
    return [0, 0, 0]

def B_sim_simple(x, y, z, B=1):
    r = np.sqrt(x**2+z**2)
    y0 = 1 - np.abs(y)
    zr = 1 - 0.5 * np.abs(z)
    #vec = np.array([-z, 0, x])
    scale = np.nan_to_num(r/y0**2 * (1-np.abs(r/y0 - 1)), nan=0.0, posinf=0.0, neginf=0.0)
    return [B * scale * -z, 0, B * scale * x]

from scipy.integrate import quad

def R(a, b, c):
    return np.matrix([
        [np.cos(a)*np.cos(b),
         np.cos(a)*np.sin(b)*np.sin(c)-np.sin(a)*np.cos(c),
         np.cos(a)*np.sin(b)*np.cos(c)+np.sin(a)*np.sin(c)],

        [np.sin(a)*np.cos(b),
         np.sin(a)*np.sin(b)*np.sin(c)+np.cos(a)*np.cos(c),
         np.sin(a)*np.sin(b)*np.cos(c)-np.cos(a)*np.sin(c)],

        [-np.sin(b), np.cos(b)*np.sin(c), np.cos(b)*np.cos(c)]
    ])

def P_unint(x, y, z, a=0, b=0, c=0, B=1):
    R_inv = R(a, b, c).getT()
    inp = np.array(R_inv.dot([x, y, z]))[0]
    B_sim_true = B_sim(*inp, B)
    B_rot = np.array(R(a, b, c).dot(B_sim_true))[0]
    return B_rot[2]

def P_int(x, y, a=0, b=0, c=0, B=1):
    lamb = lambda z: P_unint(x, y, z, a, b, c, B)
    result = quad(lamb, -10, 10, limit=500)
    return 0.5 * result[0]

print("Part 1 Complete")

# === PART 2 - CREATE BACKGROUND === #

from PIL import Image
import numpy as np
import random
import math
from scipy.stats import norm

#def plot_RMs(rms, scale=1, rm_scale=1):
#        rms_pos = np.array(list(filter(lambda val: val[2] > 0, rms.transpose()))).transpose()
#        rms_neg = np.array(list(filter(lambda val: val[2] < 0, rms.transpose()))).transpose()
#
#        plt.scatter(rms_pos[0] * scale, rms_pos[1] * scale, s=rms_pos[2]*rm_scale, color=(1, 1, 1, 0), edgecolors='red', linewidth=2)
#        plt.scatter(rms_neg[0] * scale, rms_neg[1] * scale, s=-1 * rms_neg[2]*rm_scale, color=(1, 1, 1, 0), edgecolors='blue', linewidth=2)

def create_background_grid(resolution=30, gridn=2, scale=256, bias=128, plot=False):
    # Create random RM grid
    points = resolution*gridn**2

    image = Image.open('../foreground/cmb_cropped_BW.jpg')
 
    # Summarize some details about the image
    #print(image.format)
    #print(image.size)
    #print(image.mode)
    #image.show()

    # Get the image data as an array
    width, height = image.size
    image_arr = list(map(lambda px: int((px[0]+px[1]+px[2])/3), list(image.getdata())))
    image_arr = [image_arr[i*width:(i+1)*width] for i in range(height)]
    #print(image_arr)

    # Generate simulated lists
    xy = np.array(random.sample([(x,y) for x in range(width) for y in range(height)], points))
    x = xy.transpose()[0]
    y = xy.transpose()[1]
    s = []

    # Develop 30 random numbers as POIs
    for i in range(points):
        # Set zero point to be 50% of 256, or 128
        s.append((image_arr[x[i]][y[i]] - bias)/scale)

    # Convert to single degree scale
    x = gridn * x/(width)
    y = gridn * y/(height)

    ss = np.array([x, y, s])

    #if plot:
    #    plt.imshow(image)
    #    plot_RMs(ss, width/gridn, scale)
    #    plt.show()

    return ss

print("Part 2 Complete")

#create_background_grid(int(np.rint(30*size**2)), bias=150, plot=True)

# === PART 3 - SAVE ALL MODELS === #

from astropy.table import Table
from astropy.io import ascii

def write_processed(table, file):
    table.write(file+".ecsv", format='ascii.ecsv', overwrite=True)

def read_processed(file):
    return ascii.read(file+".ecsv")

# Create background variants
profiles = [100, 100, 128, 128, 150, 150]
profs = []
for p in range(len(profiles)):
    size = 0.25 * (1 + (random.random() * (np.pi - 1)))
    x, y, r = create_background_grid(int(np.rint(30*size**2)), gridn=8, bias=profiles[p])
    t = Table([x-4, y-4, r, np.sqrt(np.abs(r)), np.zeros(len(r)), np.zeros(len(r))], names=('x','y','RM','RM_uncert', 'interpolation_raw', 'interpolation_unc'))
    inner = t[((x-4)**2+(y-4)**2)<4]
    outer = t[((x-4)**2+(y-4)**2)>4]
    print(len(inner),len(outer))
    write_processed(outer, "../data_processed/toy_model/background_models/outer_"+str(p+1))
    write_processed(inner, "../data_processed/toy_model/background_models/inner_"+str(p+1))
    profs.append(inner)

import copy

angles = np.array([0,45,90,135,180,225,270,315]) * np.pi/180
alpha = angles
beta = angles
gamma = angles
opacity = [0,0.2,0.4,0.6,0.8,1] # Remove opacity from evaluations?

# Code for file naming:
# toy_[background]_[alpha]_[beta]_[gamma]

master_list = Table(names=["background", "alpha", "beta", "gamma"])

long = len(profiles) * len(angles) ** 3

nc = 0
for a in alpha:
    for b in beta:
        for c in gamma:
            for p in range(len(profs)):
                master_list.add_row([p, a, b, c])
                inner_prof = copy.deepcopy(profs[p])
                for i in range(len(inner_prof)):
                    bigP = P_int(inner_prof[i]["x"], inner_prof[i]["y"], a, b, c)
                    inner_prof[i]["RM"] = inner_prof[i]["RM"] + bigP
                    inner_prof[i]["RM_uncert"] = np.sqrt(np.abs(inner_prof[i]["RM"]) + np.abs(bigP))
                write_processed(inner_prof, "../data_processed/toy_model/toy_hvcs/toy_"+str(p)+"_"+str(a)+"_"+str(b)+"_"+str(c))
                del inner_prof
                nc = nc + 1
                print("Creating models: "+str(int(100*nc/long))+"% (Model#"+str(nc)+"/"+str(long)+") \r", sep="", end="", flush=True)

write_processed(master_list, "../data_processed/toy_model/master_model")

print("Part 3 Completed")
print("Termination")