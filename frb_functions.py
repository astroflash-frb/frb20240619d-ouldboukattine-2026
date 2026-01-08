import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chi2
from scipy import integrate
from scipy.optimize import minimize
from scipy import special
from astropy.coordinates import SkyCoord, EarthLocation, Angle, AltAz
from astropy.units import Quantity
from tqdm import tqdm
import copy
import re
import your
import pandas as pd
import astropy.units as u


from typing import Optional, List, Dict, Union, Tuple

def chime_coord_converter(source_name, alpha_deg, delta_deg, output_format="both"):
    """
    Converts celestial coordinates (RA, Dec) from degrees to formatted strings in various styles.

    Parameters:
    -----------
    source_name : str
        The name of the (FRB) source.
    alpha_deg : float
        Right Ascension (RA) in degrees. Must be between 0 and 360.
    delta_deg : float
        Declination (Dec) in degrees. Must be between -90 and +90.
    output_format : str, optional
        The format in which to output the coordinates. Options are:
        - "key": Format suitable for .key file
        - "wb": Format suitable for Hour Angle Wb script
        - "both": Output both formats (default)
    
    Returns:
    --------
    dict
        A dictionary containing the formatted RA and Dec.
        
    Examples:
    ---------
    Convert coordinates for a source and print both key and Wb script formats:
    chime_coord_converter(source_name="R49", alpha_deg=8.265, delta_deg=+28.831, wb_fmt=True, key_fmt=True)
    """
    # Validate input ranges
    if not (0 <= alpha_deg <= 360):
        raise ValueError("Right Ascension (alpha_deg) must be between 0 and 360 degrees.")
    
    if not (-90 <= delta_deg <= 90):
        raise ValueError("Declination (delta_deg) must be between -90 and +90 degrees.")
    
    # Create SkyCoord object
    source = SkyCoord(ra=alpha_deg * u.degree, dec=delta_deg * u.degree, frame='icrs')
    
    # Convert to formatted strings
    ra_hour_fmt = source.ra.to_string(u.hour, sep=':', pad=True)
    dec_hour_fmt = source.dec.to_string(u.deg, sep=':', pad=True, alwayssign=True)

    output = {}
    output["ra"] = ra_hour_fmt
    output["dec"] = dec_hour_fmt

    if output_format in ["key", "both"]:
        print("format for in .key file:")
        print(f"source='{source_name}' ra={ra_hour_fmt} dec={dec_hour_fmt} equinox='j2000' / \n")

    if output_format in ["wb", "both"]:
        print("format for in Hour Angle Wb script")
        print(f"'{source_name}' : SkyCoord('{ra_hour_fmt}', '{dec_hour_fmt}', unit=(u.hourangle, u.deg)), \n")

    return output
    
def calc_sweep(chans, tsamp, dm, dmconst=1.0 / (2.41 * 10.0**-1.0), printres=True):
    """ Function to calculate how long the dispersive sweep is in units of
    milliseconds and in units of time bins 
                   
    Input:         
    chans: 1d numpy array containing channel frequencies in units of MHz
    dm: float, dispersion measure in units of pc cm^-3
    tsamp: float, sampling time in units of seconds
    dmconst: float, the dispersion constant in units of GHz^2 cm^3 pc^-1 ms
        default = 1.0 / (2.41 * 10.0**-1.0) = 4.149377593360995 

    Returns:       
    time_of_sweep_ms: float, the duration of the sweep in milliseconds
    bin_sweep: the duration of the sweep in units of the sampling time
    """            
                   
    # calculate sweep duration in ms
    time_of_sweep_ms = dmconst * 10**6.0 * dm * (np.min(chans)**-2.0 - np.max(chans)**-2.0)
                   
    # calculate sweep duration in units of time bins
    # first convert time_of_sweep_ms to seconds since tsamp was given in units of seconds
    bin_sweep = np.ceil( (time_of_sweep_ms / 1000.0) / tsamp).astype(np.int64)
                   
    if printres:
        print(f"The duration of the sweep is: {(time_of_sweep_ms / 1000.):.3f} seconds")
        print(f"The duration of the sweep is: {(time_of_sweep_ms):.3f} milliseconds")
        print(f"The duration of the sweep is: {(time_of_sweep_ms * 1000.):.3f} mircoseconds")
        print(f"This corresponds to {bin_sweep} time-bins.")
    
    return time_of_sweep_ms, bin_sweep
    
def get_frequncies_from_archive(archive):
    """ Docstring """

    nchan = archive.get_nchan()
    bw = archive.get_bandwidth()
    fcent = archive.get_centre_frequency()
    cw = bw / nchan
    fmin = fcent - np.abs(bw / 2.0)
    fmax = fcent + np.abs(bw / 2.0)

    chans = archive.get_frequencies()

    return fmin, fmax, fcent, bw, cw, nchan, chans


def get_channel_centers(fmin, fmax, nchan):
    bw = fmax - fmin
    foff = bw / nchan
    chans = np.arange(nchan) * foff + (fmin + foff / 2.0)
    chans = np.flip(chans) # now the first element contains the highest frequency
    return chans
        
def downsample(data, freqf, timef, verbose=False):
    """
    Downsamples a 2D array in time (x-axis) and/or frequency (y-axis)
    Removes any time samples (before downsampling) to make the downsampling
    in time work, if needed.

    Parameters
    ----------
    data : 2d numpy array
        array containing the data, usually one of the coherence products
        (AA BB CR CI) or one of the Stokes Parameters (I Q U V).
    freqf : int
        Frequency downsample factor (use 1 if you do not want
        to downsample in frequency). Numbers of channels should be divisible
        by this factor.
    timef : int
        Time downsample factor (use 1 if you do not want to downsample in time).

    Returns
    -------
    data_down : 2d numpy array
        downsampled version of the input data array.

    """

    nchan = data.shape[0]
    nsamps = data.shape[1]

    assert nchan % freqf == 0

    remainder = nsamps % timef
    if remainder > 0:
        if verbose:
            print(f"Warning: removing the last {remainder} samples of the time-axis (before downsampling")
        data = data[:, :-remainder]

    if verbose:
        print(f'Old shape (before downsampling): {data.shape}')

    # downsample in time
    data_down = data.reshape(-1, data.shape[-1]//timef, timef).mean(axis=2)

    # downsample in frequency
    data_down = data_down.reshape(-1, freqf, data_down.shape[-1]).mean(axis=1)

    if verbose:
        print(f'New shape (after downsampling): {data.shape}')

    return data_down

def dedisperse(data, dm, freqs, tsamp, dmconst=1.0 / (2.41 * 10.0**-4.0)):
    """
    Parameters
    ----------
    data : 2d numpy array
        data containing freq on the y-axis and time on the x-axis
    dm : float
        Dispersion measure in units of pc / cc.
    freqs : 1d numpy array
        Array containing the channel frequencies in units of MHz
    tsamp : float
        sampling time in units of seconds.
    dmconst: float
        dispersion constant in units of MHz^2 pc^-1 cm^3 s
        according to K Nimmo (https://ui.adsabs.harvard.edu/abs/2022NatAs...6..393N/abstract)
        this is also the value that is used in dspsr (digifil) and SFXC

    Returns
    -------
    ndata : 2d numpy
        dedispersed data

    """
    # set up freq info
    freqs = freqs.astype(np.float64)
    reffreq = np.max(freqs)
    
    # init empty array to store dedisp data in
    ndata = np.zeros(data.shape)
    
    # calculate time shifts and convert to bin shifts
    time_shift = dmconst * dm * (reffreq**-2.0 - freqs**-2.0)
    
    # round to nearest integer
    bin_shift = np.rint((time_shift / tsamp)).astype(np.int64)
    
    # checks
    assert len(bin_shift) == data.shape[0]
    assert ndata.shape == data.shape
    
    # dedisperse by rolling back the channels
    for i, bs in enumerate(bin_shift):
        ndata[i, :] = np.roll(data[i, :], bs)
        
    return ndata

def downsample_1d(arr, ds, verbose=True):
    """
    Averages a 1D numpy array using a downsampling factor ds. Throws away potential remainders
    at the end of the 1D numpy array. If arr is a pandas series, it will first convert it to a numpy array.

    Parameters
    ----------
    arr : 1D numpy array
        Array that is to be downsampled
    ds : int
        Downsampling factor.

    Returns
    -------
    arr_down : 1e numpy array
        downsampled array.

    """
    if isinstance(arr, pd.Series):
        arr = arr.to_numpy()

    remainder = arr.shape[0] % ds
    if remainder > 0:
        if verbose:
            print(f"Warning: removing the last {remainder} samples of the time-series")
            print(f'Old shape (before downsampling): {arr.shape}')
        arr = arr[:-remainder]
        if verbose:
            print(f'Old shape (before downsampling): {arr.shape}')

    arr_down = np.mean(arr.reshape(-1, ds), axis=1)
    return arr_down    

def correct_bandpass(data, scale=3, return_extra=False, sides="both"):
    """
    Parameters
    ----------
    data : 2d numpy array
        data that needs be normalised. Usually contains frequency on the y-axis and time samples on the x-axis.
    scale : int or float, optional
        Region over which to calculate the mean and std. The default is 3.
        I.e. a scale of 3 and a sides="both" means that the first 1/3 and the last 1/3 of the data will be used
        to calculate the mean and std (i.e. it assumes that the burst is in the middle)
    return_extra : Bool, optional
        If true will also return the mean and std of every channel. The default is False.
    sides : str, optional
        if "left" will use the first 1/scale part of the data to calculate the statistics.
        if "right" will use the last 1/scale part of the data to calculate the statistics.
        if "both" will use both sides
        The default is "both".

    Returns
    -------
    2d numpy array
        data that has been bandpass corrected, i.e. (data - mean) / std (on a per channel basis)
        if return_extra = True will also return two 1d numpy arrays containing the mean and std of every channel

    """

    assert sides in ["left", "right", "both"]

    nsamp = data.shape[1]
    idx = np.arange(nsamp)
    idx_off_left = idx[:nsamp//scale]
    idx_off_right = idx[-nsamp//scale:]

    if sides == "left":
        idx_off = idx_off_left
    elif sides == "right":
        idx_off = idx_off_right
    elif sides == "both":
        idx_off = np.concatenate((idx_off_left, idx_off_right))

    bp_mean = np.mean(data[:, idx_off], axis=1)
    bp_std = np.std(data[:, idx_off], axis=1)

    data_bpc = (data - bp_mean[:, None]) / bp_std[:, None]

    if return_extra:
        return data_bpc, bp_mean, bp_std

    return data_bpc

def normalize_1d(data, scale=3, sides='both'):
    """ Function that normalized an 1d array.
    Defines the 'offburst' region as the left 1/scale part of the data and
    the rightern most 1/scale part of the data. """

    assert sides in ["left", "right", "both"]

    nsamp = data.shape[0]
    idx = np.arange(nsamp)
    idx_off_left = idx[:nsamp//scale]
    idx_off_right = idx[-nsamp//scale:]

    if sides == "left":
        idx_off = idx_off_left
    elif sides == "right":
        idx_off = idx_off_right
    elif sides == "both":
        idx_off = np.concatenate((idx_off_left, idx_off_right))

    data_off_mean = np.mean(data[idx_off])
    data_off_std = np.std(data[idx_off])

    data_norm = (data - data_off_mean) / data_off_std
    return data_norm
        
def is_numeric(variable):
    # Check if the variable is an instance of int, float, or any numpy integer or float type
    return isinstance(variable, (int, float, np.integer, np.floating))
    
def mask_array(array, bad_chan_idx):
    """
    Masks entire rows in a 2D NumPy array.

    Converts the input array into a masked array and masks the rows specified 
    in `bad_chan_idx`, ignoring them in computations.

    Parameters:
    -----------
    array : numpy.ndarray
        A 2D NumPy array.
    bad_chan_idx : list of int
        Indices of rows to be masked.

    Returns:
    --------
    numpy.ma.MaskedArray
        The input array with specified rows masked.

    Example:
    --------
    >>> arr = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    >>> mask_array(arr, [1])
    masked_array(
      data=[[1, 2, 3],
            [--, --, --],
            [7, 8, 9]],
      mask=[[False, False, False],
            [ True,  True,  True],
            [False, False, False]]
    )
    """
    temp_array = np.ma.masked_array(array, mask=False)
    temp_array.mask[bad_chan_idx] = True
    return temp_array


def mask_dict(data_dict, bad_chan_idx):
    dict_temp = {}
    for k, v in data_dict.items():
        temp_array = np.ma.masked_array(v, mask=False)
        temp_array.mask[bad_chan_idx] = True
        dict_temp[k] = mask_array(v, bad_chan_idx)
    return dict_temp

def get_subband_edges(bandwidth, n_chan, subbands=8):
    
    # ab is starting array [0,1,126,127]
    subband_edge = n_chan // subbands 
    
    #if statement so to make sure then bandpass also gets smaller when downsampling
    ab = np.array([0,1,subband_edge-2,subband_edge-1])
        
    if n_chan <= 512:
        ab = np.array([0,subband_edge-1])
        
    print("bandwidth is:", bandwidth, ". Number of channels:", n_chan, ". First subband array: ", ab)
    print("Subband edges: ", subband_edge)
        
    mask = []
    
    for i in range(8):
        a = ab+(i*subband_edge)
        mask.append(a)
            
    mask_arr = np.asarray(mask).ravel()

    return mask_arr
        
def your_header_to_dict(yr_obj):
    """ To do: write docstring """

    header_dict = {"basename": yr_obj.your_header.basename,\
                   "bw": yr_obj.your_header.bw,\
                   "center_freq": yr_obj.your_header.center_freq,\
                   "dec_deg": yr_obj.your_header.dec_deg,\
                   "dtype": yr_obj.your_header.dtype,\
                   "fch1": yr_obj.your_header.fch1,\
                   "filename": yr_obj.your_header.filename,\
                   "foff": yr_obj.your_header.foff,\
                   "format": yr_obj.your_header.format,\
                   "frequency_decimation_factor": yr_obj.your_header.frequency_decimation_factor,\
                   "gb": yr_obj.your_header.gb,\
                   "gl": yr_obj.your_header.gl,\
                   "native_foff": yr_obj.your_header.native_foff,\
                   "native_nbits": yr_obj.your_header.native_nbits,\
                   "native_nchans": yr_obj.your_header.native_nchans,\
                   "native_nspectra": yr_obj.your_header.native_nspectra,\
                   "native_tsamp": yr_obj.your_header.native_tsamp,\
                   "nbits": yr_obj.your_header.nbits,\
                   "nchans": yr_obj.your_header.nchans,\
                   "npol": yr_obj.your_header.npol,\
                   "nspectra": yr_obj.your_header.nspectra,\
                   "poln_order": yr_obj.your_header.poln_order,\
                   "ra_deg": yr_obj.your_header.ra_deg,\
                   "source_name": yr_obj.your_header.source_name,\
                   "time_decimation_factor": yr_obj.your_header.time_decimation_factor,\
                   "tsamp": yr_obj.your_header.tsamp,\
                   "tstart": yr_obj.your_header.tstart,\
                   "tstart_utc": yr_obj.your_header.tstart_utc,\
                   "freqs": yr_obj.chan_freqs,\
                   "fmin": yr_obj.your_header.center_freq - np.abs(yr_obj.your_header.bw) / 2.0,\
                   "fmax": yr_obj.your_header.center_freq + np.abs(yr_obj.your_header.bw) / 2.0,\
                   }

    return header_dict
    
def get_data(yourobj, nstart=0, nsamp=None, npoln=1,\
             keys="infer"):
    """
    Parameters
    ----------
    yourobj : your object
        your object that contains metadata of a given fits or fil file.
    nstart : int
        starting sample of the data that is wanted
        default: 0
    nsamp : None or int
        if None will default to all time samples
        default: None (all samples)
    npoln : int
        Number of polarisations. Should be 1 (intensity) or 4 (full pol or coherence)
        default: 1
    keys : string or list with strings
        keys to use in the dict that is returned. If this is "infer" and if npoln = 4 will use ["AA", "BB", "CR", "CI"].
        default: "infer"

    Returns
    -------
    If npoln = 1 returns a 2D numpy array
    If npoln = 4 returns a dictonary with four keys, each pointing to their own 2D numpy array
    """

    assert npoln == 1 or npoln == 4

    if keys == "infer" and npoln == 4:
        keys = ["AA", "BB", "CR", "CI"]

    if npoln == 4:
        assert type(keys) == list and len(keys) == npoln

    if nsamp == None:
        nsamp = int(yourobj.native_nspectra)

    data = yourobj.get_data(nstart=nstart, nsamp=nsamp, npoln=npoln).T

    if npoln == 1:
        return data

    res = {}
    # store data in dict
    for i in range(npoln):
        res[keys[i]] = data[:, i, :]

    return res

def extract_numbers_pazi(input_string, map_channels=False, total_channels=None):
    """
    Extracts space-separated numbers from a pazi string enclosed in double quotes. 
    Optionally, remaps channels such that channel 0 maps to 127, 1 to 126, etc. in case you want it flipped

    Args:
        input_string (str): The input string containing space-separated numbers in double quotes.
        map_channels (bool): If True, remaps channels based on the total_channels argument (default is False).
        total_channels (int, optional): Total number of channels, required if map_channels is True.

    Returns:
        np.ndarray: A numpy array of extracted and optionally remapped integers.
    
    Example:
        input_string = "paz -z \"49 50 45 46 28 30 15 127 16 103\" -e pazi B01-NRT.fil.ar"
        result = extract_numbers_pazi(input_string, map_channels=True, total_channels=128)
        print(result)  # Output: [49 50 45 46 28 30 15 127 16 103]
    """
    if map_channels and total_channels is None:
        raise ValueError("total_channels must be provided if map_channels is True.")
    
    match = re.search(r'"([\d\s]+)"', input_string)
    if match:
        numbers_str = match.group(1)
        numbers = np.array([int(num) for num in numbers_str.split()])
        
        if map_channels:
            max_channel = total_channels - 1
            numbers = max_channel - numbers
        
        return numbers
    return np.array([])

