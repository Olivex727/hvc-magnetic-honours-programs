#!/usr/bin/env python3

import numpy as np

from astropy.io import ascii

from astropy.wcs import WCS
from astropy.io import fits

from astropy.coordinates import SkyCoord
from astropy.io.votable import parse_single_table

from foreground import foreground_remover

class file_find:

    def get_RMs():
        return ascii.read("../data_catalog/CASDA_RMs.ecsv")

    def get_HI_emission():
        hdu = fits.open("../data_catalog/hi4pi-hvc-nhi-ait.fits")
        header = hdu[0].header
        image = hdu[0].data
        wcs = WCS(header)
        return image, header
    
    def get_HVC_locations(hvc_area_range=(1, np.pi), full_hvc_range=False):
        def format_sexagesimal(coord):
            ra, dec = coord.split()[:3], coord.split()[3:]
            ra_formatted = ':'.join(ra)
            dec_formatted = f"{dec[0]}:{dec[1]}:{dec[2]}"
            return f"{ra_formatted} {dec_formatted}"

        # Read the VOTable file
        table = parse_single_table('../data_catalog/vizier_Moss2013_HVCs.vot').to_table()

        # Mask to decent sized area
        table_big_HVCs = table
        if (not full_hvc_range):
            mask = (table['Area']>hvc_area_range[0]) & (table['Area']<hvc_area_range[1])
            table_big_HVCs = table[mask]

        # Assuming the coordinates are in columns 'ra' and 'dec', adjust as needed
        ra_dec = [f"{row['RAJ2000']} {row['DEJ2000']}" for row in table_big_HVCs]
        ra_dec = [format_sexagesimal(f"{row['RAJ2000']} {row['DEJ2000']}") for row in table_big_HVCs]

        ra_dec = SkyCoord(ra_dec,frame='fk5',unit=('hour','deg')).icrs

        # Data is now in desired structure
        table.add_column(ra_dec, index=1, name="SkyCoord")

        return table

    def get_H_alpha():
        return 0
    
    def get_corrected_RMs(h_alpha, rm_grid, catalog="POSSUM"):
        return foreground_remover.subtract_foreground(h_alpha, rm_grid, catalog)

class collator:

    def collate_whole_sky(interpolation_catalog, hvc_area_range=(1, np.pi), full_hvc_range=False):
        RMs = file_find.get_RMs()
        H_alpha = file_find.get_H_alpha()
        return RMs, file_find.get_HI_emission(), file_find.get_HVC_locations(hvc_area_range, full_hvc_range), H_alpha, file_find.get_interpolation(H_alpha, RMs)