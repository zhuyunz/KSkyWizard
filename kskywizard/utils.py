from pypeit.par import pypeitpar
from pypeit.core import telluric
from pypeit.spectrographs.util import load_spectrograph
import zap
import os
import re
import numpy as np
from scipy.interpolate import interp1d
from astropy.io import fits
import math
import ref_index
import pkg_resources
try:
    from ConfigParser import ConfigParser
except ImportError:
    from configparser import ConfigParser

def read_setup_cfg(section, key):
    config = ConfigParser()

    # Get the directory where this script resides
    setup_cfg_path = os.path.join(os.path.dirname(__file__), '..', 'setup.cfg')

    # Read the setup.cfg file
    config.read(setup_cfg_path)

    # Access the value in custom_section
    my_value = config.get(section, key)

    return my_value

telgridfile = read_setup_cfg('processing', 'telgridfile')
try:
    hdu = fits.open(telgridfile)
except:
    raise ValueError(f"Telgrid file: '{telgridfile}', in 'setup.cfg' does not exist.")

def telluric_correct(infile_path: str,  star_ra: float, star_dec: float, 
                     spectrograph = 'keck_kcrm'):
    """
    Telluric correct the spectrum of a given standard star.
    Original written by Milan Sharma Mandigo-Stoba for DBSP_DRP at https://github.com/finagle29/DBSP_DRP/blob/main/dbsp_drp/telluric.py;
    Adatpted by Zhuyun Zhuang for the KCWI data

    Args:
        coadd (str): Coadd filename.
        star_ra (float), star_dec (float): the RA and Dec of the STD
        spectrograph (str): PypeIt name of spectrograph.
    """
    spectrograph = load_spectrograph(spectrograph)
    par = spectrograph.default_pypeit_par()

    #only used for standard star
    par['telluric']['objmodel'] = 'star'
    par['telluric']['star_ra'] = star_ra
    par['telluric']['star_dec'] = star_dec
    # if par['telluric']['telgridfile'] is None:
    #     if par['sensfunc']['IR']['telgridfile'] is not None:
    #         par['telluric']['telgridfile'] = par['sensfunc']['IR']['telgridfile']
    par['telluric']['telgridfile'] = telgridfile
    par['telluric']['teltype'] = 'grid'

    # Parse the output filename
    outfile = re.sub('.fits', '_tellcorr.fits', infile_path)
    modelfile = re.sub('.fits', '_tellmodel.fits', infile_path)

    try:
        TelStar = telluric.star_telluric(infile_path, par['telluric']['telgridfile'],
                                        modelfile, outfile,
                                        star_type=par['telluric']['star_type'],
                                        star_mag=par['telluric']['star_mag'],
                                        star_ra=par['telluric']['star_ra'],
                                        star_dec=par['telluric']['star_dec'],
                                        func=par['telluric']['func'],
                                        model=par['telluric']['model'],
                                        polyorder=par['telluric']['polyorder'],
                                        only_orders=par['telluric']['only_orders'],
                                        teltype=par['telluric']['teltype'], tell_npca=par['telluric']['tell_npca'],
                                        mask_hydrogen_lines=par['sensfunc']['mask_hydrogen_lines'],
                                        mask_helium_lines=par['sensfunc']['mask_helium_lines'],
                                        hydrogen_mask_wid=par['sensfunc']['hydrogen_mask_wid'],
                                        delta_coeff_bounds=par['telluric']['delta_coeff_bounds'],
                                        minmax_coeff_bounds=par['telluric']['minmax_coeff_bounds'],
                                        pix_shift_bounds=par['telluric']['pix_shift_bounds'],
                                        maxiter=par['telluric']['maxiter'],
                                        popsize=par['telluric']['popsize'],
                                        tol=par['telluric']['tol'])
    except ValueError:
        print(f"[ERROR] Telluric correction of {os.path.base(infile_path)} FAILED!")


def kcwi_correct_extin(img0, hdr0):
    """Atmospheric extinction correction from official KCWI DRP"""
    img = img0.copy()
    hdr = hdr0.copy()

    # get airmass
    air = hdr['AIRMASS']

    # read extinction data
    full_path = pkg_resources.resource_filename(__name__, 'data/extin/snfext.fits')
    if os.path.exists(full_path):
        hdul = fits.open(full_path)
        exwl = hdul[1].data['LAMBDA']
        exma = hdul[1].data['EXT']
        # get object wavelengths
        sz = img.shape
        dw = hdr['CD3_3']
        w0 = hdr['CRVAL3']
        owls = np.arange(sz[0]) * dw + w0
        # linear interpolation
        exint = interp1d(exwl, exma, kind='cubic', bounds_error=False,
                         fill_value='extrapolate')
        # resample extinction curve
        oexma = exint(owls)
        # convert to flux ratio
        flxr = 10.**(oexma * air * 0.4)
        if len(sz) == 3:
            # apply to cube
            for ix in range(sz[2]):
                for iy in range(sz[1]):
                    img[:, iy, ix] *= flxr
        else:
            # apply to vector
            img *= flxr

        flrmn = np.nanmean(flxr)
        hdr['HISTORY'] = 'kcwi_correct_extin'
        hdr['EXTCOR'] = (True, 'extinction corrected?')
        hdr['AVEXCOR'] = (flrmn, 'average extin. correction (flux ratio)')
        # if logger:
        #     logger.info("Extinction corrected")
        # else:
        print("Extinction corrected")
            
        return img, hdr
    else:
        print("Extinction data file (%s) not found!" % full_path)

def scale_extinct_sky(hdr_sky, hdr_sci):
    """Atmospheric extinction correction from official KCWI DRP"""

    # get airmass
    air_sky = hdr_sky['AIRMASS']
    air_sci = hdr_sci['AIRMASS']
    # read extinction data
    full_path = pkg_resources.resource_filename(__name__, 'data/extin/snfext.fits')
    if os.path.exists(full_path):
        hdul = fits.open(full_path)
        exwl = hdul[1].data['LAMBDA']
        exma = hdul[1].data['EXT']
        # get object wavelengths
        dw = hdr_sky['CD3_3']
        w0 = hdr_sky['CRVAL3']
        owls = np.arange(hdr_sky['NAXIS3']) * dw + w0
        # linear interpolation
        exint = interp1d(exwl, exma, kind='cubic', bounds_error=False,
                         fill_value='extrapolate')
        # resample extinction curve
        oexma = exint(owls)
        # convert to flux ratio
        flxr_sky = 10.**(oexma * air_sky * 0.4)
        flxr_sci = 10.**(oexma * air_sci * 0.4)
            
        return flxr_sky / flxr_sci
    
    else:
        print("Extinction data file (%s) not found!" % full_path)


def reassign_skyseg(x_segments, widths, new_segment, new_width, min_span = 100):
    """
    Reassign the cfwidth to skyseg required by MUSE.  Primarily written by ChatGPT.
    Using skyseg defined by MUSE directly would cause weird fluctuations if a strong emission line of the science object
    falls in the region where the sky lines are sparse. This function is to replace the cfwidth of the
    regions w/ strong emission lines to a narrower value

    Parameters:
    x_segments (arr, with shape of n): An array of segment grid points.
    widths: (arr, with shape of n-1): The cfwidth of each segment.
    new_segment (arr): The new segment as [a1, a2].
    new_width (float): The weight for the new segment.
    min_span (int): the minimial length of each segment

    Returns:
    list of tuple: The updated list of segments with their weights.
    """
    a1, a2 = new_segment
    a1 = int(a1)
    a2 = int(a2)

    # Add new segment to the list of segments and sort
    x_segments.extend([a1, a2])
    x_segments = sorted(set(x_segments))

    # Initialize new widths list
    updated_widths = []

    # Assign widths to the new segments
    for i, (start, end) in enumerate(zip(x_segments[:-1], x_segments[1:])):
        if a1 <= start < a2:
            updated_widths.append(new_width)
        else:
            width = widths[min(i, len(widths) - 1)]  # Ensure valid index access
            updated_widths.append(width)

    # Merge adjacent segments and ensure each spans at least min_span
    final_segments = []
    final_widths = []

    i = 0
    while i < len(x_segments) - 1:
        start, end = x_segments[i], x_segments[i + 1]
        width = updated_widths[i]

        # Merge adjacent segments
        while i < len(x_segments) - 1 and end - start < min_span:
            next_start, next_end = x_segments[i + 1], x_segments[i + 2]
            if next_start - end < min_span:
                end = next_end
                width = new_width
                i += 1
            else:
                break

        final_segments.append((start, end))
        final_widths.append(width)
        i += 1

    # Convert final_segments to 1D array
    final_x_segments = []
    for start, end in final_segments:
        if not final_x_segments or final_x_segments[-1] != start:
            final_x_segments.append(start)
        final_x_segments.append(end)

    return final_x_segments, final_widths


def atm_disper(w0, w1, airmass, temperature=10.0, pressure_pa=61100.0,
               humidity=50.0, co2=400.0):
    """

    Calculate atmospheric dispersion at w1 relative to w0 [Copied from KCWI_DRP, put it here to avoid re-running the DRP...]

    Args:
        w0 (float): reference wavelength (Angstroms)
        w1 (float): offset wavelength (Angstroms)
        airmass (float): unitless airmass
        temperature (float): atmospheric temperature (C)
        pressure_pa (float): atmospheric pressure (Pa)
        humidity (float): relative humidity (%)
        co2 (float): Carbon-Dioxide (mu-mole/mole)

    """

    # Calculate
    z = math.acos(1.0/airmass)

    n0 = ref_index.ciddor(wave=w0/10., t=temperature, p=pressure_pa,
                          rh=humidity, co2=co2)
    n1 = ref_index.ciddor(wave=w1/10., t=temperature, p=pressure_pa,
                          rh=humidity, co2=co2)

    return 206265.0 * (n0 - n1) * math.tan(z)

def collapse_header(hdr):
    """
    Quick wrapper to collapse a 3-D header into a 2-D one.
    Copied from KCWIKit

    Parameters
    ----------
    hdr: header

    Returns
    -------
    hdr_img: collapsed header

    """

    hdr_img=hdr.copy()
    hdr_img['NAXIS']=2
    del hdr_img['NAXIS3']
    del hdr_img['CD3_3']
    del hdr_img['CTYPE3']
    del hdr_img['CUNIT3']
    del hdr_img['CNAME3']
    del hdr_img['CRVAL3']
    del hdr_img['CRPIX3']

    return hdr_img

def check_dir(dir):
    if not os.path.isdir(dir):
        os.mkdir(dir)

    return
