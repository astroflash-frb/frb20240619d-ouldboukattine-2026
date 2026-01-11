import your
import frb_functions
import argparse
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import matplotlib

from funcs import autocorr, acf_per_subband, fit_subband_acfs_direct, fit_fullband_acf_direct
     

# -----------------------------
# Command-line argument parser
# -----------------------------
parser = argparse.ArgumentParser(description="Process FRB burst dynamic spectrum.")
parser.add_argument("file", type=str, help="Name or path of the burst file.")
parser.add_argument("--downsamp_f", type=int, default=1,
                    help="Frequency downsampling factor (default: 1)")
parser.add_argument("--downsamp_t", type=int, default=1,
                    help="Time downsampling factor (default: 1)")
parser.add_argument("--maxlag", type=float, default=10,
                    help="Max lag in MHz for ACF (default: 10MHz)")
parser.add_argument(
    "--noise",
    action="store_true",
    help="Use this flag if you want to compute the off-burst ACF"
)

args = parser.parse_args()
burst = args.file
filename=burst
downsamp_f = args.downsamp_f
downsamp_t = args.downsamp_t
maxlag = args.maxlag
noise = args.noise
burst = your.Your(burst)

#header information

header_dict = frb_functions.your_header_to_dict(burst)

tsamp = header_dict["tsamp"]
nbins = header_dict["nspectra"]
fsamp = header_dict["foff"] * -1
freqs = header_dict["freqs"] # in Mhz

tmin = 0.0
tmax = nbins * burst.tsamp # seconds
fmin = header_dict["fmin"]
fmax = header_dict["fmax"]
extent_b = [tmin, tmax, fmin, fmax]

#2D numpy array
d = frb_functions.get_data(burst, nstart=0, nsamp=nbins, npoln=1)

#bandpass correction
sides_bandpass = 'both'
scale_bandpass = 3
stokes_i = frb_functions.correct_bandpass(d, sides=sides_bandpass, scale=scale_bandpass)
#stokes_i=stokes_i[:,100*32:250*32]

#downsample
stokes_i = frb_functions.downsample(stokes_i, freqf=downsamp_f, timef=downsamp_t)

#time profile
prof = np.nanmean(stokes_i,axis=0)
stokes_i_normprof = frb_functions.normalize_1d(prof, sides=sides_bandpass, scale=scale_bandpass)

plt.plot(stokes_i_normprof)
plt.show()

# Plot dynamic spectrum
"""
fig, ax_dyn = plt.subplots(figsize=(8, 6))
ax_dyn.imshow(stokes_i, aspect='auto', interpolation='none')#,extent=extent)
ax_dyn.set_xlabel("Time [s]")
ax_dyn.set_ylabel("Frequency [MHz]")
plt.show()
"""

if 'B05' in filename:
    #spec = (np.nanmean(stokes_i[:,15803:15850], axis=1) + np.nanmean(stokes_i[:,np.argmax(prof):np.argmax(prof)+2],axis=1)) /2.
    spec = stokes_i[:,np.argmax(prof)]
    if noise:
        #offspec = np.nanmean(stokes_i[:,:49],axis=1)
        offspec = stokes_i[:,49]
elif 'B07' in filename:
    #spec = (np.nanmean(stokes_i[:,11804:11806],axis=1) + stokes_i[:,np.argmax(prof)] + stokes_i[:,11800] + np.nanmean(stokes_i[:,11854:11896],axis=1)) / 4.
    spec = stokes_i[:,np.argmax(prof)]
    if noise:
        offspec = stokes_i[:,50]
        #offspec = np.nanmean(stokes_i[:,:46],axis=1)
elif 'B11_p55081_19_128us' in filename:
    #spec = (np.nanmean(stokes_i[:,5292:5299],axis=1) + np.nanmean(stokes_i[:,5309:5314],axis=1)) /2.
    spec = np.nanmean(stokes_i[:,5310:5314],axis=1)
    if noise:
        offspec = np.nanmean(stokes_i[:,500:504],axis=1)
elif 'B06_p55049_25_256us_2KHz_StokesI_DM465095_scin.cor_Wb' in filename:
    #spec = (np.nanmean(stokes_i[:,np.argmax(prof):np.argmax(prof)+2],axis=1) + stokes_i[:,1179]) / 2.
    spec = stokes_i[:,np.argmax(prof)]
    if noise:
        offspec = stokes_i[:,200]
elif 'p55025_28_64us_7KHz_StokesI_DM464857.cor_Wb.fil' in filename:
    spec = (np.nanmean(stokes_i[:,5142:5255],axis=1) + np.nanmean(stokes_i[:,5295:5316],axis=1) + np.nanmean(stokes_i[:,5560:5567],axis=1)) / 3.
elif 'B01_p55025_28_128us_4KHz_StokesI_DM464857_scin.cor_Wb' in filename:
    spec = stokes_i[:,np.argmax(prof)]
    if noise:
        offspec = stokes_i[:,200]

elif 'B01_p55025_28_8ms_61Hz_StokesI_DM464857_scin.cor_Wb' in filename:
    #spec = np.nanmean(stokes_i[:,41:43],axis=1)
    spec = stokes_i[:,np.argmax(prof)]
    if noise:
        offspec = stokes_i[:,20]

elif 'B10' in filename:
    spec = stokes_i[:,np.argmax(prof)]#np.nanmean(stokes_i[:,298:304],axis=1)
    if noise:
        offspec = stokes_i[:,6]
        #offspec = np.nanmean(stokes_i[:,:6],axis=1)
else:
    spec=stokes_i[:,np.argmax(prof)]

if 'B02' or 'B06_p55049_25_4ms_122Hz_StokesI_DM465095_scin_Wb' or 'B11_p55081_19_4ms_122Hz_StokesI_DM465266_scin.cor_Wb' in filename:
    if noise:
        offspec = stokes_i[:,5]

def downsample_freqs(freqs, downsample):
    # Ensure freqs is a NumPy array
    freqs = np.array(freqs)

    # Trim to nearest multiple of `downsample` if needed
    trimmed_len = len(freqs) - (len(freqs) % downsample)
    freqs_trimmed = freqs[:trimmed_len]

    # Reshape and average along the downsample axis
    freqs_ds = freqs_trimmed.reshape(-1, downsample).mean(axis=1)

    return freqs_ds

if downsamp_f !=1:
    fsamp*=downsamp_f
    downsamp_freqs = downsample_freqs(freqs, downsamp_f)
    freqs=downsamp_freqs

if noise:
    plt.plot(freqs,offspec,color='k')

plt.plot(freqs,spec,color='k',alpha=0.3)
plt.show()


#now we have our spectrum lets ACF
maxlag_bin = int(maxlag//fsamp)
offspec_mean = np.nanmean(stokes_i[:,0])
acf=autocorr(spec, v=None,zerolag=False,maxlag=maxlag_bin,offspec_mean=offspec_mean,freq=None)
if noise:
    acf_noise = autocorr(offspec, v=None,zerolag=False,maxlag=maxlag_bin,offspec_mean=offspec_mean,freq=None)
height=acf[0]
lags=np.arange(len(acf))+1

acf=acf[1:]
lags=lags[1:]
acf=np.concatenate((acf[::-1],acf))
lags=np.concatenate((-1*lags[::-1],lags))*fsamp
if noise:
    acf_noise=acf_noise[1:]
    acf_noise=np.concatenate((acf_noise[::-1],acf_noise))

results = fit_fullband_acf_direct(acf, lags, fsamp, 
                lagrange_for_fit=5, 
                force_model=None,
                xlim=maxlag)

if noise:
    try:
        results_noise = fit_fullband_acf_direct(acf_noise, lags, fsamp,
                lagrange_for_fit=5,
                force_model=None,
                xlim=maxlag)
    except:
        results_noise=None
# Save
bname=filename[8:11]
if noise:
    np.savez('acf_results_%s_wnoise_new.npz'%bname, acf=acf, lags=lags, acf_noise=acf_noise,results=results, results_noise=results_noise)
else:
    np.savez('acf_results_%s_new.npz'%bname, acf=acf, lags=lags, results=results)

#plt.plot(lags,acf,drawstyle='steps-mid',color='k',linewidth=0.5)
#plt.show()

spec[0]=0
spec=np.ma.masked_where(spec==0,spec)

numsubs=16

acfs,fcents,lags, subsn, submask,spec_lens=acf_per_subband(spec,freqs,np.linspace(0,len(freqs)-1,len(freqs)),num_subbands=numsubs,plot_fit=False,savefig=None,maxlag=maxlag,offspec=stokes_i[:,0],snsubband=True)
cmap = matplotlib.cm.get_cmap('plasma')
for i in range(len(fcents)):
        rgba = cmap(i/len(fcents))
        plt.plot(lags[len(fcents)-i-1],acfs[len(fcents)-i-1]+(1.0*i),drawstyle='steps-mid',color=rgba,linewidth=2,alpha=0.7,label='%.2f MHz'%fcents[len(fcents)-i-1])
        plt.xlabel('Freq Lag [MHz]')
        plt.ylabel('ACF power [offset + 1]')

plt.show()

results = fit_subband_acfs_direct(
    acfs=acfs,
    fcents=fcents,
    lags=lags,
    f_res = fsamp,
    lagrange_for_fits=[0.025]*numsubs,  # compute out to 2MHz
    offset=1.0,                            # y-axis plotting offset
    xlim=0.02,#maxlag,                              # fit only within ± maglag MHz
    force_model=None                      # try both single and double Lorentzian
)
