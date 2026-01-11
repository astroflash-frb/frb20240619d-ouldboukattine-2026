import numpy as np
import matplotlib.pyplot as plt
import os
from funcs import lorentz_w_c, doublelorentz_w_c

# === Parameters ===
burst_ids = ['B01','B02', 'B05', 'B06', 'B07', 'B10', 'B11']
largerscale=True
vals = ['d_l','s','s','s','s','s','s']

#burst_ids = ['B01','B02','B06','B10', 'B11']
#vals=['d_l','s','s','s','s']

xlim =1.0  # or whatever value you typically use
ncols = 7
nrows = 2#int(np.ceil(len(burst_ids) / ncols))

# === Plotting ===
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(5 * ncols, 4 * nrows), squeeze=False)

for i, bid in enumerate(burst_ids):
    # Load data
    if largerscale is True:
        file=f'acf_results_{bid}_wnoise.npz'
        wnoise=True
    else:
        try:
            file = f'acf_results_{bid}_highres.npz'
            wnoise=False
            if not os.path.exists(file):
                print(f"File not found: {file}")
                print(something)
        except:
            file = f'acf_results_{bid}_highres_wnoise.npz'
            wnoise=True
            if not os.path.exists(file):
                print(f"File not found: {file}")
                continue

    data = np.load(file, allow_pickle=True)
    acf = data['acf']
    lag = data['lags']
    if wnoise is True:
        acf_noise = data['acf_noise']
        results_noise = data['results_noise'].item()

    results = data['results'].item()  # unpack dictionary

    #print(bid)
    #print(results)

    # === Helper function to plot ===
    def plot_acf(ax, lag, acf, results, label_prefix):
        ax.plot(lag, acf, drawstyle='steps-mid', color='gray', lw=2, label='Data')

        # 1-Lorentzian fit
        try:
            p1 = results['1_lorenz']
            model_vals1 = lorentz_w_c(lag, gamma1=p1['scint'], m1=p1['mod'], c=p1['c'])
            ax.plot(lag, model_vals1, 'r--', label='1-Lorentzian')

            # 2-Lorentzian fit
            p2 = results['2_lorenz']
            model_vals2 = doublelorentz_w_c(
                lag,
                gamma1=p2['scint_1'], m1=p2['mod_1'],
                gamma2=p2['scint_2'], m2=p2['mod_2'],
                c=p2['c']
        )
            ax.plot(lag, model_vals2, 'k--', label='2-Lorentzian')
        
            # Extract appropriate scint and mod based on vals

            if vals[i] is 's':
                min_scint=p1['scint']
                scint_err=p1['scint_uncert']
                min_mod = p1['mod']
                mod_err = p1['mod_uncert']
            elif vals[i] is 'd_l':
                if p2['scint_1'] < p2['scint_2']:
                    min_scint = p2['scint_2']
                    scint_err=p2['scint_uncert_2']
                    min_mod = p2['mod_2']
                    mod_err = p2['mod_uncert_2']
                else:
                    min_scint = p2['scint_1']
                    scint_err=p2['scint_uncert_1']
                    min_mod = p2['mod_1']
                    mod_err = p2['mod_uncert_1']
            elif vals[i] is 'd_s':
                if p2['scint_1'] < p2['scint_2']:
                    min_scint = p2['scint_1']
                    scint_err=p2['scint_uncert_1']
                    min_mod = p2['mod_1']
                    mod_err = p2['mod_uncert_1']
                else:
                    min_scint = p2['scint_2']
                    scint_err=p2['scint_uncert_2']
                    min_mod = p2['mod_2']
                    mod_err = p2['mod_uncert_2']

            min_scint *= 1000  # to kHz
            scint_err*=1000
            freq_res = np.abs(lag[1] - lag[0]) * 1000  # to kHz
            
            #additional uncertainty
            bw=128*1000
            N = 1 + 0.2 * (bw / min_scint)
            add_un = min_scint / (2 * np.sqrt(N))

            print("bw: %s+-%s  kHz"%(min_scint,np.sqrt(scint_err**2+add_un**2)))
            print("m: %s+-%s"%(min_mod,mod_err))#np.sqrt(mod_err**2+add_un**2)))
            
            ax.text(
                0.05, 0.95,
                f"$\\nu_d$ = {min_scint:.1f} kHz\n$m$ = {min_mod:.2f}\n$\\Delta\\nu$ = {freq_res:.2f} kHz",
                transform=ax.transAxes,
                va='top', ha='left', fontsize=10,
                bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='gray', alpha=0.8)
            )
        except:
            print(' ')

        ax.set_xlim(-xlim, xlim)
        ax.set_xlabel("Frequency lag [MHz]")
        ax.set_ylabel("ACF")
        ax.set_title(f"{bid} ({label_prefix})")
        #ax.legend()

    """
    ax = axes[i // ncols, i % ncols]
    ax.plot(lag, acf, drawstyle='steps-mid', color='gray', lw=2, label='Data')

    # 1-Lorentzian
    p1 = results['1_lorenz']
    model_vals1 = lorentz_w_c(lag, gamma1=p1['scint'], m1=p1['mod'], c=p1['c'])
    ax.plot(lag, model_vals1, 'r--', label='1-Lorentzian')

    # 2-Lorentzian
    p2 = results['2_lorenz']
    model_vals2 = doublelorentz_w_c(
        lag,
        gamma1=p2['scint_1'],
        m1=p2['mod_1'],
        gamma2=p2['scint_2'],
        m2=p2['mod_2'],
        c=p2['c']
    )
    ax.plot(lag, model_vals2, 'k--', label='2-Lorentzian')

    # Determine smallest scint and corresponding mod from 2 lorentz fit
    if bid == 'B01' or bid == 'B11' or 'highres' in file:
        if p2['scint_1'] < p2['scint_2']:
            min_scint = p2['scint_1']
            min_mod = p2['mod_1']
        else:
            min_scint = p2['scint_2']
            min_mod = p2['mod_2']
    else:
        #use single lorentzian vals for legend
        min_scint = p1['scint']
        min_mod = p1['mod']

    min_scint*=1000
    # Frequency resolution
    freq_res = np.abs(lag[1] - lag[0])
    freq_res*=1000
    # Add text to panel
    ax.text(
        0.05, 0.95,
        f"$\\nu_d$ = {min_scint:.1f} kHz\n$m$ = {min_mod:.2f}\n$\\Delta\\nu$ = {freq_res:.2f} kHz",
        transform=ax.transAxes,
        va='top', ha='left',
        fontsize=10,
        bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='gray', alpha=0.8)
    )

    ax.set_xlim(-xlim, xlim)
    ax.set_xlabel("Frequency lag [MHz]")
    ax.set_ylabel("ACF")
    ax.set_title(bid)
    ax.legend()
    """

    # === Top row: original ACF ===
    ax_top = axes[0, i]
    plot_acf(ax_top, lag, acf, results, "on burst")

    # === Bottom row: noise-subtracted ACF ===
    ax_bot = axes[1, i]
    if wnoise is True:
        plot_acf(ax_bot, lag, acf_noise, results_noise, "off burst")

    if nrows==3:
        ax_ex = axes[2,i]
        if wnoise is True:
            acfmax=np.max(acf)
            acf_corr = acf/acfmax - acf_noise/np.max(acf_noise)
            acf_corr*=acfmax
            plot_acf(ax_ex, lag, acf_corr, None, "corrected")

fig.tight_layout()
plt.show()
