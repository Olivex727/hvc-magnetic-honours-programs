#!/usr/bin/env python3

from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs import WCS
import numpy as np

import astropy_healpix as healpix

import matplotlib.pyplot as plt

class image_transform:
    def get_pixel(wcs, coords):
        car = 0
        if "galactic" in wcs.world_axis_physical_types[0]:
            coords = coords.galactic
            car = np.array([[coords.l.deg],[coords.b.deg]]).T
        else:
            coords = coords.icrs
            car = np.array([[coords.ra.deg],[coords.dec.deg]]).T

        return WCS.wcs_world2pix(wcs,car,1)[0]
    
    # Image must be in FITS format
    def get_flux_at_point(image, coords, healpix=False):
        if healpix:
            return 0
        else:
            wcs = WCS(image.header)
            pix = image_transform.get_pixel(wcs, coords)
            return image.data[int(pix[1])-1][int(pix[0])-1]
    
    def healpy_to_fits(healpy):
        return 0