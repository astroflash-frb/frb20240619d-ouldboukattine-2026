import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from lmfit import minimize, Parameters, fit_report, Model, Minimizer, report_fit
from tqdm import tqdm 

#necessary fcns
def shift(v, i, nchan):
        """
        function v by a shift i
        nchan is the number of frequency channels (to account for negative lag)
        """
        n = len(v)
        r = np.zeros(3*n)
        i+=nchan-1 #to account for negative lag
        i = int(i)
        r[i:i+n] = v
        return r

def autocorr(spec, v=None,zerolag=False,maxlag=None,offspec_mean=None,freq=None):
    """
    x is the 1D array you want to autocorrelate
    v is the array of 1s and 0s representing a mask where 1 is no mask, and 0 is mask
    zerolag = True will keep the zero lag noise spike, otherwise it won't compute the zero lag
    maxlag = None will compute the ACF for the entire length of x
    maxlag = bin_number will compute the ACF for lags up to x[bin_number]

    """
    nchan=len(spec)
    if v is None:
        v = np.ones_like(spec)


    x = np.copy(spec)
    xmean=np.nanmean(x[v!=0])

    if offspec_mean is None:
        denom = xmean**2
    else:
        denom = (xmean - offspec_mean)**2


    x[v!=0] -= xmean#x[v!=0].mean()
    if maxlag==None:
        ACF = np.zeros_like(x)
    else:
        ACF = np.zeros_like(x)[:int(maxlag)]

    for i in tqdm(range(len(ACF))):
        if zerolag == False:
                if i>1:
                        m = shift(v,0,nchan)*shift(v,i,nchan)
                        ACF[i-1] = np.nansum(shift(x,0,nchan)*shift(x, i,nchan)*m) / (np.sum(m)*denom)
        else:
                m = shift(v,0,nchan)*shift(v,i,nchan)
                ACF[i] = np.nansum(shift(x,0,nchan)*shift(x, i,nchan)*m) / (np.sum(m)*denom)

    return ACF

def acf_per_subband(spec,freqs,freqids,num_subbands=2,savefig='./acf_per_freq.pdf',plot_fit=True,maxlag=None,snsubband=False,offspec=None):

    plt.close()

    sub_cent=[]
    sub_scint=[]

    sub_len = len(spec)//num_subbands
    acfs=[]
    lags=[]
    fcents=[]
    sub_sn=[]
    sub_mask=[]
    spec_lens=[]
    mask = np.abs(np.array(spec.mask, dtype='bool')-1)
    tot=np.sum(spec.data*mask)
    for sub in range(num_subbands):
        if snsubband is False:
            beg = sub*sub_len
            end = (sub+1)*sub_len
            if end>(len(spec)-1):
                end=-1
            subtot=np.sum(spec.data[beg:end])

        else:
            if sub==0:
                beg = 0
            else:
                beg = end

            i=beg-1
            subtot=0
            while subtot < (tot/float(num_subbands)):
                i+=1
                subtot+=(spec.data[i]*mask[i])

            end = i

        sub_sn.append(subtot)
        sub_mask.append(np.sum(spec.mask[beg:end]))
        if end!=-1:
            spec_lens.append(end-beg)
        else:
            spec_lens.append(len(spec)-beg)

        print("beg,end",beg,end)
        if 5 >= maxlag:
            lagrange=maxlag
        else:
            lagrange=5

        if offspec is not None:
            acf = acf_scint_plot(spec[beg:end],freqids[beg:end],freqs[beg:end],[0,0],lagrange_for_fit=lagrange,diagnostic_plots=False,maxlag=maxlag,offspec_mean=np.nanmean(offspec[beg:end]))
        else:
            acf = acf_scint_plot(spec[beg:end],freqids[beg:end],freqs[beg:end],[0,0],lagrange_for_fit=lagrange,diagnostic_plots=False,maxlag=maxlag)

        acfs.append(acf[0])
        lags.append(acf[1])
        cmap = matplotlib.cm.get_cmap('plasma')
        rgba = cmap(sub/num_subbands)

        plt.plot(acf[1],acf[0]+(1*sub),drawstyle='steps-mid',color=rgba,linewidth=1,alpha=1,label='%.2f MHz'%(freqs[beg]+((freqs[end]-freqs[beg])/2)))
        fcents.append(freqs[beg]+((freqs[end]-freqs[beg])/2))
        if plot_fit==True:
            try:
                plt.plot(acf[1],lorentz(acf[1],acf[2].params['gamma'],acf[2].params['m'],acf[2].params['c']) +(1*sub),color='k',linewidth=0.5)
                sub_cent.append((freqs[beg]+((freqs[end]-freqs[beg])/2)))
                sub_scint.append(np.abs(acf[2].params['gamma']))
            except:
                sub_cent.append((freqs[beg]+((freqs[end]-freqs[beg])/2)))
                sub_scint.append(0)
    if savefig!=None:
        plt.xlim(-maxlag,maxlag)
        plt.ylim(-1,1+(sub))
        plt.xlabel('Freq lag [MHz]')
        plt.legend(loc='upper left')
        plt.savefig(savefig,format='pdf')
    else:
        plt.close('all')

    if plot_fit==True and savefig!=None:
        plt.close()
        plt.scatter(sub_cent,sub_scint,marker='x',color='k')
        plt.plot(freqs,sub_scint[-1]*(freqs/sub_cent[-1])**4,color='r')
        plt.xlabel('Freq [MHz]')
        plt.ylabel('Scint bw [MHz]')
        plt.savefig(savefig[:-4]+'_scintbw.pdf',format='pdf')

    # Ensure all ACFs (and corresponding lags) have the same length
    acf_lengths = [len(acf) for acf in acfs]
    min_len = min(acf_lengths)
    f_res = np.abs(freqs[1]-freqs[0])

    if len(set(acf_lengths)) > 1:
        print(f"Warning: Center-trimming ACFs to minimum common length {min_len//2*f_res} MHz")
        acfs_trimmed = []
        lags_trimmed = []
        for acf, lag in zip(acfs, lags):
            diff = len(acf) - min_len
            start = diff // 2
            end = start + min_len
            acfs_trimmed.append(acf[start:end])
            lags_trimmed.append(lag[start:end])
        acfs = acfs_trimmed
        lags = lags_trimmed

    return acfs,fcents,lags, sub_sn, sub_mask, spec_lens

def acf_scint_plot(ds,freq_ids,freqs,time_range,lagrange_for_fit=10.,diagnostic_plots=True, maxlag=None, offspec_mean=None):
    """
    ds is either 2D array [freq,time] or 1d spectrum
    freq_ids and freqs are the frequency channel numbers and central frequency in MHz mapping ds
    time_range is [begin_bin,end_bin] within which the burst spectrum is computed. Only needed if ds is 2d.
    lagrange_for_fit is the lag range in MHz used to define the lag range out to which the final lorentzian will be fit
    diagnostic_plots will produce plots while running
    maxlag is the maximum lag in MHz to compute the ACF out to
    """

    #figure out the frequency resolution
    num_chan_diff = int(np.abs(freq_ids[1]-freq_ids[0]))
    f_res = np.abs(freqs[1]-freqs[0])/float(num_chan_diff)
    print("Frequency resolution is %.5f MHz"%f_res)

    #make the spectrum
    if ds.ndim == 1:
        spec = ds
    else:
        spec=np.mean(ds[:,time_range[0]:time_range[1]],axis=1)

    try:
        if ds.ndim == 1:
            mask =  np.abs(np.array(ds.mask,dtype=int)-1)
        else:
            mask=np.abs(np.array(ds.mask[:,0],dtype=int)-1)
    except:
        print('masking where the array = 0')
        ds = np.ma.masked_where(ds==0,ds)
        if ds.ndim == 1:
            mask =  np.abs(np.array(ds.mask,dtype=int)-1)
        else:
            mask=np.abs(np.array(ds.mask[:,0],dtype=int)-1)

    #ACF
    if maxlag ==None:
        maxlag_bin=None
    else:
        maxlag_bin=int(maxlag/f_res)

    acf=autocorr(spec, v=mask,zerolag=False,maxlag=maxlag_bin,offspec_mean=offspec_mean,freq=None)
    height=acf[0]
    lags=np.arange(len(acf))+1

    acf=acf[1:]
    lags=lags[1:]
    acf=np.concatenate((acf[::-1],acf))
    lags=np.concatenate((-1*lags[::-1],lags))*f_res

    #plotting
    if diagnostic_plots==True:
        plt.plot(lags,acf,drawstyle='steps-mid',color='k',linewidth=0.5)
        plt.show()

        plt.plot(lags,acf,drawstyle='steps-mid',color='k',linewidth=0.5,label="%.5f MHz"%f_res)
    #fit lorentzian to measure scintillation bandwidth
    try:
        gmodel = Model(lorentz_w_c)
        acf_for_fit = acf[int(len(acf)/2.)-int(lagrange_for_fit/f_res):int(len(acf)/2.)+int(lagrange_for_fit/f_res)]
        lags_for_fit = lags[int(len(acf)/2.)-int(lagrange_for_fit/f_res):int(len(acf)/2.)+int(lagrange_for_fit/f_res)]
        result = gmodel.fit(acf_for_fit, x=lags_for_fit, gamma=0.001, m=1, c=0)
        if diagnostic_plots == True:
            plt.plot(lags,lorentz_w_c(lags,result.params['gamma'],result.params['m'],result.params['c']),color='orange',label='scint bw = %.2f MHz'%result.params['gamma'].value)
            plt.xlim(-np.abs(result.params['gamma'].value)*10,np.abs(result.params['gamma'].value)*10)
            plt.ylim(-0.2,height+0.05)
            plt.legend()
            plt.show()

        return acf, lags, result
    except:
        print("Could not fit a Lorentzian")
        if diagnostic_plots==True:
            plt.legend()
            plt.show()

        return acf, lags

def fit_subband_acfs_direct(acfs, fcents, lags, f_res,
                             lagrange_for_fits=None,
                             offset=1.0,
                             xlim=0.5,
                             force_model=None, ignore_subband_idx=None):
    """
    Fit subband ACFs with 1- and 2-Lorentzian models and compare.

    Parameters
    ----------
    acfs : list of np.ndarray
        ACFs per subband.
    fcents : list or array
        Central frequencies per subband.
    lags : list of np.ndarray
        Frequency lags per subband.
    f_res: float 
        Frequency resolution in MHz
    lagrange_for_fits : list of float, optional
        Lag fit ranges (in MHz) per subband. Defaults to 0.3 for all.
    offset : float, optional
        Vertical offset between ACF plots.
    xlim : float, optional
        X-axis limit for lag plots.
    force_model : int, optional
        Force use of 1 or 2 Lorentzians. Default is None (try both).

    Returns
    -------
    results_dict : dict
        Dictionary of fitting results.
    """
    numsubs = len(acfs)
    cmap = matplotlib.cm.get_cmap('plasma')

    if lagrange_for_fits is None:
        lagrange_for_fits = [0.3] * numsubs

    results_dict = {'1_lorenz': {}, '2_lorenz': {}, 'f_cents': [], 'acfs_offset': [], 'lm_fitting_objects': []}
    for key in ['1_lorenz', '2_lorenz']:
        results_dict[key] = {}

    fig = plt.figure(figsize=(5, 7))

    # initalize outputs for the models
    for key in list(results_dict.keys()):
        if '_lorenz' in key:
            results_dict[key] = {"sub_scint_1": [], "sub_scint_uncert_1": [], "mods1": [], "mods1_uncert": [], "add_un1": [], 'c1': [], 'c1_uncert': []}
            
            

            if int(key[0]) == 2:
                results_dict[key] = {"sub_scint_1": [], "sub_scint_uncert_1": [], "mods1": [], "add_un1": [],  "add_un2": [],
                                      "mods1_uncert": [], "sub_scint_2": [], "sub_scint_uncert_2": [], 
                                      "mods2": [], "mods2_uncert": []}

    for i in range(numsubs):
        print(f'i : {i} fcent: {fcents[i]} total_offset {float(offset)*i}')
        lag = lags[i]
        acf = acfs[i]
        rgba = cmap(i / numsubs)
        lagrange_for_fit = lagrange_for_fits[i]

        results_dict["acfs_offset"].append(acf + offset * i)
        results_dict["f_cents"].append(fcents[i])


        acf_fit = acf[int(len(acf) / 2.) - int(lagrange_for_fit / f_res):int(len(acf) / 2.) + int(lagrange_for_fit / f_res)]
        lag_fit = lag[int(len(acf) / 2.) - int(lagrange_for_fit / f_res):int(len(acf) / 2.) + int(lagrange_for_fit / f_res)]
        acf_fit_errors = calculate_acf_errors(acf_fit)

        plt.plot(lag, acf + offset * i, drawstyle='steps-mid', color=rgba, linewidth=2, alpha=0.7, label=f'{fcents[i]:.2f} MHz')

        fit_results = {}
        for numlorentz, model_func in zip(
            [1, 2], [lorentz_withc_min, doublelorentz_withc_min]
        ):
            if force_model is not None and numlorentz != force_model:
                continue

            params = Parameters()
            #params.add('gamma1', value=np.float(lagrange_for_fit), min=1e-5, max=100)
            params.add('gamma1', value=0.05, min=1e-5, max=100)
            params.add('m1', value=1, min=0, max=100)
            if numlorentz == 2:
                #params.add('gamma2', value=np.float(lagrange_for_fit)/10., min=1e-5, max=100)
                params.add('gamma2', value=0.001, min=1e-5, max=100)
                params.add('m2', value=1, min=0, max=100)
            params.add('c', value=0, min=-100, max=100)

            try:
                fit_min = Minimizer(model_func, params, fcn_args=(lag_fit, acf_fit, np.sqrt(acf_fit_errors)))
                result = fit_min.minimize()
                fit_results[numlorentz] = result
            except Exception as e:
                print(f"Fitting failed for subband {i}, model {numlorentz}: {e}")

        results_dict['lm_fitting_objects'].append(fit_results)

        # Store results for 1 Lorentzian
        result_1 = fit_results[1]
        if result_1:
            label = None
            if i == 0:
                label = 'Single'
            plt.plot(lag, lorentz_w_c(lag, result_1.params['gamma1'], result_1.params['m1'], result_1.params['c']) + (float(offset) * i), color='red', linestyle='-.', linewidth=1, label=label)
            results_dict['1_lorenz']["sub_scint_1"].append(np.abs(result_1.params['gamma1'].value))
            results_dict['1_lorenz']["sub_scint_uncert_1"].append(abs(result_1.params['gamma1'].stderr))
            results_dict['1_lorenz']["mods1"].append(np.abs(result_1.params['m1'].value))
            results_dict['1_lorenz']["mods1_uncert"].append(abs(result_1.params['m1'].stderr))
            results_dict['1_lorenz']['c1'].append(result_1.params['c'].value)
            results_dict['1_lorenz']['c1_uncert'].append(result_1.params['c'].stderr)


            scints = np.array([np.abs(result_1.params['gamma1'].value)])
            modinds = np.array([np.abs(result_1.params['m1'].value)])
            errs = np.array([abs(result_1.params['gamma1'].stderr)])
            moderrs = np.array([abs(result_1.params['m1'].stderr)])
            inds = np.argsort(scints)


            # Compute additional uncertainty from low number of scintles
            sub_scint1 = scints[inds][0]
            good_chans = np.array([len(x) for x in acfs])#np.array(spec_lens) - np.array(submask)
            #might need to change this above line if we include an RFI mask 

            N = 1 + 0.2 * ((np.flip(good_chans) * f_res) / sub_scint1)
            add_un = sub_scint1 / (2 * np.sqrt(N))

            results_dict['1_lorenz']["add_un1"].append(add_un)


        else:
            results_dict['1_lorenz']["sub_scint_1"].append(0)
            results_dict['1_lorenz']["sub_scint_uncert_1"].append(0)
            results_dict['1_lorenz']["mods1"].append(0)
            results_dict['1_lorenz']["mods1_uncert"].append(0)
            results_dict['1_lorenz']["add_un1"].append(0)
            results_dict['1_lorenz']['c1'].append(0)
            results_dict['1_lorenz']['c1_uncert'].append(0)


        # Store results for 2 Lorentzians
        result_2 = fit_results[2]
        if result_2:
            label = None
            if i == 0:
                label = 'Double'
            plt.plot(lag, doublelorentz_w_c(lag, result_2.params['gamma1'], result_2.params['m1'], result_2.params['gamma2'], result_2.params['m2'], result_2.params['c']) + (float(offset) * i), color='k', linestyle='--', linewidth=1, label=label)
            scints = np.array([np.abs(result_2.params['gamma1'].value), np.abs(result_2.params['gamma2'].value)])
            modinds = np.array([np.abs(result_2.params['m1'].value), np.abs(result_2.params['m2'].value)])

            gamma1_stderr = result_2.params['gamma1'].stderr
            gamma2_stderr = result_2.params['gamma2'].stderr
            mod1_stderr = result_2.params['m1'].stderr
            mod2_stderr = result_2.params['m2'].stderr

            # Check if stderr is None and replace it with np.nan
            errs = np.array([abs(gamma1_stderr) if gamma1_stderr is not None else np.nan,
                             abs(gamma2_stderr) if gamma2_stderr is not None else np.nan])
            # errs = np.array([abs(result_2.params['gamma1'].stderr), abs(result_2.params['gamma2'].stderr)])
            moderrs = np.array([abs(mod1_stderr) if mod1_stderr is not None else np.nan,
                             abs(mod2_stderr) if mod2_stderr is not None else np.nan])

            inds = np.argsort(scints)
            results_dict['2_lorenz']["sub_scint_1"].append(scints[inds][0])
            results_dict['2_lorenz']["sub_scint_2"].append(scints[inds][1])

            results_dict['2_lorenz']["sub_scint_uncert_1"].append(errs[inds][0])
            results_dict['2_lorenz']["sub_scint_uncert_2"].append(errs[inds][1])


            results_dict['2_lorenz']["mods1"].append(modinds[inds][0])
            results_dict['2_lorenz']["mods1_uncert"].append(moderrs[inds][0])
            results_dict['2_lorenz']["mods2"].append(modinds[inds][1])
            results_dict['2_lorenz']["mods2_uncert"].append(moderrs[inds][1])

            # Compute additional uncertainty from low number of scintles
            sub_scint2 = scints[inds][1]
            sub_scint1 = scints[inds][0]
            good_chans = np.array([len(x) for x in acfs])#np.array(spec_lens) - np.array(submask)
            #might need to change this above line if we include an RFI mask

            N2 = 1 + 0.2 * ((np.flip(good_chans) * f_res) / sub_scint2)
            add_un2 = sub_scint2 / (2 * np.sqrt(N2))


            N = 1 + 0.2 * ((np.flip(good_chans) * f_res) / sub_scint1)
            add_un = sub_scint1 / (2 * np.sqrt(N))

            results_dict['2_lorenz']["add_un1"].append(add_un)
            results_dict['2_lorenz']["add_un2"].append(add_un2)

        else:
            results_dict['2_lorenz']["sub_scint_1"].append(0)
            results_dict['2_lorenz']["sub_scint_2"].append(0)
            results_dict['2_lorenz']["sub_scint_uncert_1"].append(0)
            results_dict['2_lorenz']["sub_scint_uncert_2"].append(0)
            results_dict['2_lorenz']["mods1"].append(0)
            results_dict['2_lorenz']["mods1_uncert"].append(0)
            results_dict['2_lorenz']["mods2"].append(0)
            results_dict['2_lorenz']["mods2_uncert"].append(0)
            results_dict['2_lorenz']["add_un1"].append(0)
            results_dict['2_lorenz']["add_un2"].append(0)


        #plot
        if 1 in fit_results:
            model_vals = lorentz_w_c(lag, **{k: fit_results[1].params[k].value for k in ['gamma1', 'm1', 'c']})
            plt.plot(lag, model_vals + offset * i, color='red', linestyle='-.', linewidth=1, label='Single' if i == 0 else None)

        if 2 in fit_results:
            model_vals = doublelorentz_w_c(lag, **{k: fit_results[2].params[k].value for k in ['gamma1', 'm1', 'gamma2', 'm2', 'c']})
            plt.plot(lag, model_vals + offset * i, color='k', linestyle='--', linewidth=1, label='Double' if i == 0 else None)

    plt.yticks([acf[0] + offset * i for i, acf in enumerate(acfs)], [f'{f:.1f}' for f in fcents])
    plt.xlim(-xlim, xlim)
    plt.xlabel('Freq lag [MHz]')
    plt.ylabel('Subband freq [MHz]')
    fig.tight_layout(pad=1)
    plt.show()

    best_model = bic_comparison(results_dict, fcents, './')

    # Now we fit for powerlaw fitting 
    power_law_fits_results_list, y1000_list, y1000_list_err = perform_powerlaw_fit('./', results_dict, best_model, ignore_subband_idx=ignore_subband_idx, force_model=force_model)

    #plot the modulation indices
    mi, mi_err = plot_mod_inds('./', results_dict, best_model)

    return results_dict


def bic_comparison(results_dict, fcents_subs, outdir, save=False):
    """
    Plots the Bayesian Information Criterion (BIC) for different lmfit models across subbands.
    Determines the best model per subband and identifies the overall best model.
    
    Parameters:
        results_dict (dict): A dictionary containing 'lm_fitting_objects', a list of subband fit results.
        fcents_subs (list): A list of central frequencies corresponding to each subband (MHz).
        outdir (str): Path to outdir to store files
        save (bool, optional): If True (default), saves the plot.
        
    returns: int, the best model IE 2 if the best model shows 2 lorentzians
    """
    lmfit_objects = results_dict['lm_fitting_objects']
    num_subbands = len(lmfit_objects)
    cmap = matplotlib.cm.get_cmap('plasma')
    
    model_labels = ["1 Lorentzian", "2 Lorentzians"]
    model_indices = [1, 2]
    best_model_counts = {1: 0, 2: 0}
    
    plt.figure(figsize=(6, 5))
    
    for i, (subband, fcent) in enumerate(zip(lmfit_objects, fcents_subs)):
        rgba = cmap(i / num_subbands)
        if 1 in subband and 2 in subband:  # Ensure both models exist
            bic1 = subband[1].bic
            bic2 = subband[2].bic
            best_model = 1 if bic1 < bic2 else 2
            best_model_counts[best_model] += 1
            plt.scatter(1, bic1, color=rgba)
            plt.scatter(2, bic2, color=rgba)
            plt.text(best_model, min(bic1, bic2), f'✓', ha='center', fontsize=12, fontweight='bold', color='black')
    
    overall_best_model = max(best_model_counts, key=best_model_counts.get)
    plt.axvline(x=overall_best_model, linestyle="--", color="k", alpha=0.7, label=f"Overall Best Model: {model_labels[overall_best_model - 1]}")
    
    plt.xticks(model_indices, model_labels)
    plt.xlabel("Model")
    plt.ylabel("BIC Value")
    plt.title(f"BIC Comparison of Models Across Subbands")
    sm = plt.cm.ScalarMappable(cmap=cmap)
    sm.set_array(fcents_subs)
    cbar = plt.colorbar(sm, label="Subband Central Frequency (MHz)")
    plt.legend()
     
    plt.show()
    
    return overall_best_model

def doublelorentz_w_c(x,gamma1,m1,gamma2,m2,c):
        return m1**2 / (1+(x/gamma1)**2) + m2**2 / (1+(x/gamma2)**2) +c

def lorentz_w_c(x,gamma1,m1,c):
        return m1**2 / (1+(x/gamma1)**2) + c

def triplelorentz(x,gamma1,m1,gamma2,m2,gamma3,m3):
        return m1**2 / (1+(x/gamma1)**2) + m2**2 / (1+(x/gamma2)**2) + m3**2 / (1+(x/gamma3)**2)

def lorentz(x,gamma1,m1):
        return m1**2 / (1+(x/gamma1)**2)

def lorentz_withc_min(params,x,y,err):
        gamma1 = params['gamma1'].value
        m1 = params['m1'].value
        c = params['c'].value

        modulo= m1**2 / (1+(x/gamma1)**2) +c
        return (modulo-y)/err

def doublelorentz_withc_min(params,x,y,err):
        gamma1 = params['gamma1'].value
        m1 = params['m1'].value
        gamma2 = params['gamma2'].value
        m2 = params['m2'].value
        c = params['c'].value

        modulo= m1**2 / (1+(x/gamma1)**2) + m2**2 / (1+(x/gamma2)**2) +c
        return (modulo-y)/err


def triplelorentz_min(params,x,y,err):
        gamma1 = params['gamma1'].value
        m1 = params['m1'].value
        gamma2 = params['gamma2'].value
        m2 = params['m2'].value
        gamma3 = params['gamma3'].value
        m3 = params['m3'].value

        modulo= m1**2 / (1+(x/gamma1)**2) + m2**2 / (1+(x/gamma2)**2) + m3**2 / (1+(x/gamma3)**2)
        return (modulo-y)/err

def calculate_acf_errors(acf):
    """
    Compute the error estimates for an autocorrelation function (ACF) using the standard method
    based on cumulative variance.

    The errors are computed using the variance formula:
        var_f[1:] = (1 + 2 * cumulative sum of squared ACF values)
    and the final errors are obtained as the square root of this variance. The function
    returns a symmetric array of errors by mirroring the computed values.

    Parameters
    ----------
    acf : array-like
        The autocorrelation function (ACF) values.

    Returns
    -------
    f_errors_full : numpy.ndarray
        The symmetric error estimates for the ACF, computed based on the variance of
        the second half of the ACF.

    """
    # Take the second half of the ACF
    acf_half = acf[len(acf) // 2:]

    # Initialize the variance for each element
    var_f = np.ones(len(acf_half)) / len(acf_half)

    # Compute the cumulative sum of squared ACF values and update var_f
    var_f[1:] *= 1 + 2 * np.cumsum(acf_half[:-1] ** 2)

    # Compute the errors as the square root of var_f
    f_errors = np.sqrt(var_f)

    # Concatenate the reversed errors with the original errors
    f_errors_full = np.concatenate((f_errors[::-1], f_errors))

    return f_errors_full


def perform_powerlaw_fit(outdir, result_dict, best_model, force_model=None,  save_plot=True, print_report=True, ignore_subband_idx= None):
    """
    Perform power-law fitting to scintillation data, including optional ignoring of subbands and plotting the results.

    This function will loop over the number of Lorentzians indicated by 'best_model', and peform a powerlaw fiting for each
    Lorentzian.


    Args:
        outdir (str): The output directory where the plots and results will be saved.
        result_dict (dict): A dictionary containing the Lorentzian model parameters and frequencies for the fitting process.
        best_model (int): The best model to use for fitting (indicates how many lorentzians to consider: 1 or 2 is supported).
        save_plot (bool, optional): Whether to save the plot of the power-law fit. Default is True.
        print_report (bool, optional): Whether to print the fitting report. Default is True.
        ignore_subband_idx (list, optional): A list of indices to exclude from the fitting process. Default is None.
        force_model (int, optional): If a value is given force fitting to use that number of Lorentzians, 1 or 2 supported, default is None.

    Returns:
        tuple: A tuple containing three elements:
            - power_law_fits_lm_results_list (list): A list of the lmfit parameter results of the power-law fitting for each lorentzian result.
            - y600_list (list): A list of the scintillation measurements at 600 MHz for each lorentzian.
            - y600_list_err (list): A list of the errors on the scintillation measurements at 600 MHz for each lorentzian.

    Example:
        >>> perform_powerlaw_fit(event_id='event123', fftsize=2048, downfreq=1, outdir='./output/', result_dict=my_result_dict, best_model=1)

    """

    if force_model:
        best_model=force_model

    f_cents = result_dict['f_cents']



    for key in (result_dict.keys()):
        if str(best_model) in key:
            best_lorentzian_params  = result_dict[key]

    params = Parameters()
    params.add('n', value=2, min=0, max=20)
    params.add('c', value=0.1)

    # we want to loop over and fit for n powerlaws
    # our model is i+1 not i

    power_law_fits_lm_results_list = [] # init
    y1000_list = [] # scintillation at 1Ghz
    y1000_list_err = [] #error on the scintillation measurement

    for i in range(int(best_model)):
        f_cents = result_dict['f_cents']
        i+=1
        sub_scint = best_lorentzian_params[f'sub_scint_{i}']
        sub_scint_uncert = best_lorentzian_params[f'sub_scint_uncert_{i}']
        add_un = best_lorentzian_params[f'add_un{i}'][i-1]


        # Apply filtering if ignore_subband_idx is provided
                # Create a mask to exclude ignored subbands
        if ignore_subband_idx is not None:

            # Ensure f_cents is a NumPy array
            f_cents = np.array(f_cents)
            sub_scint = np.array(sub_scint)
            sub_scint_uncert = np.array(sub_scint_uncert)
            add_un = np.asarray(add_un)

            ignore_subband_idx = np.array(ignore_subband_idx)
            all_indices = np.arange(len(f_cents))
            mask = ~np.isin(all_indices, ignore_subband_idx)  # Keep everything *not* in ignore_subband_idx

            print(mask)

            f_cents = f_cents[mask]
            sub_scint = sub_scint[mask]
            sub_scint_uncert = sub_scint_uncert[mask]
            add_un = add_un[mask]


        print(sub_scint)
        print(sub_scint_uncert)
        print(f_cents)
        print(add_un)

        fit_min = Minimizer(scint_freq_relation_min, params, fcn_args=(f_cents, sub_scint, np.sqrt(np.array(sub_scint_uncert) ** 2 + add_un ** 2)))
        result_scint_i = fit_min.minimize()

        # plot
        y1000, y1000_err = plot_power_law_fit('/data1/kenzie/omar_meerkat_repeater/', i, f_cents, result_scint_i, sub_scint, np.asarray(sub_scint_uncert), np.asarray(add_un))
        if print_report ==True:
            print(report_fit(result_scint_i))

        power_law_fits_lm_results_list.append(result_scint_i)
        y1000_list.append(y1000)
        y1000_list_err.append(y1000_err)

    return power_law_fits_lm_results_list, y1000_list, y1000_list_err

def plot_power_law_fit(outdir, model_name, f_cents, lmfit_params, sub_scint, sub_scint_uncert, add_un, save_plot=True):
    """
    Plots the power-law fit to scintillation bandwidth vs central frequency and overlays the fitted model with error bars.

    Parameters:
    outdir (str): The directory where the plot will be saved.
    model_name (str): Name of the model used for fitting (e.g., 'Lorentzian 1').
    f_cents (array-like): Array of central frequencies (in MHz).
    lmfit_params (lmfit.Parameters): Parameters obtained from the lmfit model fitting of the power law (includes 'c' and 'n').
    sub_scint (array-like): Measured scintillation values per subband (in kHz).
    sub_scint_uncert (array-like): Uncertainties in the measured scintillation values per subband (in kHz).
    add_un (array-like): Additional uncertainties to be added in quadrature (in kHz).
    save_plot (bool, optional): Whether to save the plot. Defaults to True.

    Returns:
    tuple: A tuple containing:
        - y_600 (float): The fitted value of the scintillation bandwidth at 600 MHz (in kHz).
        - sigma_y_600 (float): The uncertainty in the scintillation bandwidth at 600 MHz (in kHz).

    """
    # Compute total error
    total_uncert = np.sqrt(sub_scint_uncert**2 + add_un**2)

    # Extract parameters and their uncertainties
    c = lmfit_params.params['c'].value
    n = lmfit_params.params['n'].value
    sigma_c = lmfit_params.params['c'].stderr
    sigma_n = lmfit_params.params['n'].stderr


    print(f'fcent order is {f_cents}')


    # Compute fitted model
    x_axis_model = np.linspace(np.min(f_cents), np.max(f_cents), 100)
    fitted_y = c * ( x_axis_model ** n)

    # Plot data with error bars
    plt.errorbar(f_cents, np.array(sub_scint)*1000, yerr=np.array(total_uncert)*1000, fmt='o', label='Data', capsize=5)

    # Compute and report value at x = 600
    y_1000 = c * (1271 ** n) *1000


    # use residuals of our scintillation bw measurements on the power law fit to determine the uncertainties
    fit_vals = c * (np.array(f_cents) ** n)
    resids = np.array(sub_scint) - fit_vals
    sigma_y_1000 = np.nanstd(resids) * 1000

    # Overlay fitted model
    plt.plot(x_axis_model, fitted_y*1000, 'r--',
             label=f'Fit: $y = c \cdot x^n$\nAt 1.271GHz: $y = {y_1000:.3f} \pm {sigma_y_1000:.3f}$ kHz\n$n = {n:.2f} \pm {sigma_n:.2f}$')


    # Labels and legend
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Scintillation Bandwidth (kHz)')
    plt.legend()

    if save_plot == True:
        plt.savefig(f'{outdir}/{model_name}_Lorentzian_Powerlaw.pdf')
    plt.show()

    return y_1000, sigma_y_1000


def plot_mod_inds(outdir, result_dict, best_model, save_plot=False):
    """
    Plot and optionally save the Lorentzian modulation indices for each subband in a given model.

    This function retrieves the modulation indices (MI) and their uncertainties from the provided Lorentzian model parameters,
    calculates the average and standard deviation of the MI for each subband, and generates error bar plots for the MI as a
    function of frequency. The function also saves the plots to the specified output directory if requested.

    Args:
        outdir (str): The directory where the generated plots will be saved.
        result_dict (dict): A dictionary containing the Lorentzian model parameters, including the modulation indices and their uncertainties.
        best_model (int): The number of Lorentzians to process.
        save_plot (bool, optional): Whether to save the generated plots. Default is False.

    Returns:
        tuple: A tuple containing two lists:
            - mi_avg (list): The average modulation index for each subband.
            - mi_std (list): The standard deviation of the modulation index for each subband.

    """
    f_cents = result_dict['f_cents']

    for key in (result_dict.keys()):
        if str(best_model) in key:
            best_lorentzian_params  = result_dict[key]

    mi_avg = []
    mi_std = []
    alpha=[1,0.3]
    for i in range(best_model):
        i+=1
        mod_ind = best_lorentzian_params[f'mods{i}']
        mi_avg.append(np.nanmean(mod_ind))
        mi_std.append(np.nanstd(mod_ind))

        mod_ind_uncert = best_lorentzian_params[f'mods{i}_uncert']


        plt.errorbar(f_cents,mod_ind,yerr=mod_ind_uncert,marker='x', color='k', alpha=alpha[i-1], label='%s Lorentz'%(i))
    plt.title("%s Lorentzian Modulation Indices"%(i))
    plt.xlabel('Frequency [MHz]')
    plt.ylabel('Modulation index')
    plt.ylim(0,3)
    plt.axhline(1,linestyle='--',color='k',alpha=0.5)
    plt.legend()
    if save_plot == True:
        plt.savefig(f'{outdir}/{i}_Lorentzian_modind.png')

    plt.show()

    print(mi_avg,mi_std)
    return mi_avg, mi_std


def scint_freq_relation_min(params,x,y,err):
    c = params['c'].value
    n = params['n'].value

    modulo=c*(x**n)
    return (modulo-y)/err

def scint_freq_relation(v,c,n):
    return c*(v**n)


def fit_fullband_acf_direct(acf, lag, f_res, lagrange_for_fit=0.3,
                             force_model=None, xlim=0.5):
    """
    Fit a full-band ACF with 1- and 2-Lorentzian models and compare.

    Parameters
    ----------
    acf : np.ndarray
        Full-band autocorrelation function.
    lag : np.ndarray
        Frequency lags.
    f_res : float
        Frequency resolution in MHz.
    lagrange_for_fit : float, optional
        Lag fit range (in MHz). Default is 0.3.
    force_model : int, optional
        Force use of 1 or 2 Lorentzians. Default is None (try both).
    xlim : float, optional
        X-axis limit for plot.

    Returns
    -------
    results_dict : dict
        Dictionary of fit results.
    """
    results_dict = {
        '1_lorenz': {}, '2_lorenz': {}, 'acf': acf, 'lag': lag,
        'lm_fitting_object': None
    }

    acf_fit = acf[int(len(acf) / 2.) - int(lagrange_for_fit / f_res):
                  int(len(acf) / 2.) + int(lagrange_for_fit / f_res)]
    lag_fit = lag[int(len(acf) / 2.) - int(lagrange_for_fit / f_res):
                  int(len(acf) / 2.) + int(lagrange_for_fit / f_res)]
    acf_fit_errors = calculate_acf_errors(acf_fit)

    fit_results = {}
    for numlorentz, model_func in zip(
        [1, 2], [lorentz_withc_min, doublelorentz_withc_min]
    ):
        if force_model is not None and numlorentz != force_model:
            continue

        params = Parameters()
        #params.add('gamma1', value=np.float(lagrange_for_fit), min=1e-8, max=100)
        
        params.add('gamma1', value=0.05, min=1e-8, max=100)
        params.add('m1', value=1, min=0, max=100)
        if numlorentz == 2:
            #params.add('gamma2', value=np.float(lagrange_for_fit)/10., min=1e-8, max=100)
            params.add('gamma2', value=0.001, min=1e-8, max=100)
            params.add('m2', value=1, min=0, max=100)
        params.add('c', value=0, min=-100, max=100)

        try:
            fit_min = Minimizer(model_func, params, fcn_args=(lag_fit, acf_fit, np.sqrt(acf_fit_errors)))
            result = fit_min.minimize()
            fit_results[numlorentz] = result
            print(str(numlorentz)+" Lorentzians ")
            print(report_fit(result))
        except Exception as e:
            print(f"Fitting failed for model {numlorentz}: {e}")

    results_dict['lm_fitting_object'] = fit_results

    # Plot
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(lag, acf, drawstyle='steps-mid', color='gray', lw=2, label='Data')

    if 1 in fit_results:
        result = fit_results[1]
        model_vals = lorentz_w_c(lag, **{k: result.params[k].value for k in ['gamma1', 'm1', 'c']})
        ax.plot(lag, model_vals, 'r--', label='1-Lorentzian')
        results_dict['1_lorenz'] = {
            "scint": np.abs(result.params['gamma1'].value),
            "scint_uncert": abs(result.params['gamma1'].stderr),
            "mod": np.abs(result.params['m1'].value),
            "mod_uncert": abs(result.params['m1'].stderr),
            "c": result.params['c'].value,
            "c_uncert": result.params['c'].stderr,
        }

    if 2 in fit_results:
        result = fit_results[2]
        model_vals = doublelorentz_w_c(lag, **{k: result.params[k].value for k in ['gamma1', 'm1', 'gamma2', 'm2', 'c']})
        ax.plot(lag, model_vals, 'k--', label='2-Lorentzian')

        scints = np.array([np.abs(result.params['gamma1'].value), np.abs(result.params['gamma2'].value)])
        modinds = np.array([np.abs(result.params['m1'].value), np.abs(result.params['m2'].value)])
        try:
            errs = np.array([abs(result.params['gamma1'].stderr), abs(result.params['gamma2'].stderr)])
            moderrs = np.array([abs(result.params['m1'].stderr), abs(result.params['m2'].stderr)])
        except:
            errs = np.array([np.nan, np.nan])
            moderrs = np.array([np.nan, np.nan])
        inds = np.argsort(scints)

        results_dict['2_lorenz'] = {
            "scint_1": scints[inds][0],
            "scint_2": scints[inds][1],
            "scint_uncert_1": errs[inds][0],
            "scint_uncert_2": errs[inds][1],
            "mod_1": modinds[inds][0],
            "mod_uncert_1": moderrs[inds][0],
            "mod_2": modinds[inds][1],
            "mod_uncert_2": moderrs[inds][1],
            "c": result.params['c'].value,
            "c_uncert": result.params['c'].stderr,
        }

    ax.set_xlim(-xlim, xlim)
    ax.set_xlabel("Frequency lag [MHz]")
    ax.set_ylabel("ACF")
    ax.legend()
    fig.tight_layout()
    plt.show()

    return results_dict

