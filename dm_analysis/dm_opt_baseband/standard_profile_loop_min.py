#09-Sept-2024
#OuldBoukattine
import numpy as np
import pandas as pd 
import glob, argparse, pathlib
import os 

#Loading in the burst and extracting the data is done via YOUR and scripts originally made by M.Snelders
import your
import frb_functions

from tqdm import tqdm
from matplotlib import transforms


def options():
    parser = argparse.ArgumentParser(description='Script that will read a filterbank, standardize the dynamic spectrum and return the SN of a specific component in an array.', )
    general = parser.add_argument_group('Arguments')
    general.add_argument('-f', '--filterbank', type=str, required=True,
                         help='Required, Filterbank file as input')    
    general.add_argument('-dm', '--dm', type=str, required=True,
                         help='DM-value for this loop, will be appended to the csv file')  
    general.add_argument('-t', '--time', type=float, required=True,
                         help='Time in seconds where the burst is located')  
    return parser.parse_args()

def csv_maker():
    #Simple script to only create an empty csv file in my format

    #Retrieve where the .py is excecuted and save the .csv there
    script_dir = os.path.dirname(os.path.realpath(__file__))

    #Name is unique for each run
    name = script_dir + "/peak_sn_notch.csv"

    if not os.path.exists(name):
        df=pd.DataFrame(columns=["fil_name", "DM", "tsamp", "fsamp", "sn_peak", "sn_peak_idx"])
        df.to_csv(name, sep=",", header=True, index=True,)
        print(df)
    
    return name

def csv_filler(csv_name, csv_list):
    "Function to open the .csv file and append the numbers to the df and save as csv file again"

    df = pd.read_csv(csv_name, sep=",", header=0, index_col=0, na_values='NA')
    
    #In the new version of Pandas, the append method is changed to _append. 
    #You can simply use _append instead of append, i.e., df._append(df2).
    df = df._append(pd.Series(csv_list, index=df.columns[:len(csv_list)]), ignore_index=True)
        
    print(df)
    df.to_csv(csv_name, sep=",", header=True, index=True)
    print("list was appended to the csv file")
    
    
def fil_looper(fil, burst_time, dm):
    
    #Loading in the data using YOUR
    print("loading in the filterbank")
    burst = your.Your(fil)

    header_dict = frb_functions.your_header_to_dict(burst)
    tsamp = header_dict["tsamp"]
    nbins = header_dict["nspectra"]
    fsamp = header_dict["foff"] * -1
    
    #Retrieving the 2D-array
    stokes_i = frb_functions.get_data(burst,nstart=0, nsamp=nbins, npoln=1)
    
    #Correcting for the bandpass
    #Only using the right side, since the actual burst could be in the first 1/3 of the data
    stokes_i_norm = frb_functions.correct_bandpass(stokes_i, scale=3, sides="right")
    
    #Getting the time-series and normalizing, again using the right side, want to be consistent
    prof = np.mean(stokes_i_norm,axis=0)
    stokes_i_normprof = frb_functions.normalize_1d(prof,sides="right", scale=3)

    #Retrieving the bin where the peak of the burst is expected, this is not bin-precise, but could be if your inner nerd desires it.
    burst_bin = int(burst_time / tsamp)
    
    #Now take a 1ms windows around the expected peak to retrieve the sn_max
    window_bins_total = 0.001 / tsamp 
       
    burst_window_left = burst_bin - np.floor((window_bins_total/2)).astype(int)
    burst_window_right = burst_bin + np.ceil((window_bins_total/2)).astype(int)
    burst_window = stokes_i_normprof[burst_window_left:burst_window_right] 
    
    #Retrieving the min, around the nothch, S/N and the idx of the peak S/N
    sn_max = np.min(burst_window)
    sn_max_idx = np.argwhere(stokes_i_normprof == sn_max)[0][0]
    
    csv_list = [fil, dm, tsamp, fsamp, sn_max, sn_max_idx]
    print(fil, tsamp, fsamp, sn_max, sn_max_idx)

    return csv_list

def main(args):

    csv_name = csv_maker()
    print(args.filterbank)
    csv_list = fil_looper(fil = args.filterbank, burst_time = args.time, dm = args.dm)
    csv_filler(csv_name, csv_list)

if __name__ == "__main__":
    args=options()
    main(args)
