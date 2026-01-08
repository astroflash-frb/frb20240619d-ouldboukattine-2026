# ReadMe

Github repo for paper: "A HyperFlash and ÉCLAT view of the local environment and energetics of the repeating FRB 20240619D"
This repository contains the scripts and databases used to generate the plots and results presented in [Ould-Boukattine+2026](https://arxiv.org/abs/2509.16374).
The data is hosted at Zenodo can be retrieved [here](https://doi.org/10.5281/zenodo.18184163).

Software used:

The pipeline written to process and search the baseband/VDIF data can be found at [frb-baseband](https://github.com/pharaofranz/frb-baseband). software used to transfer radio data can be found at [jive5ab](https://github.com/jive-vlbi/jive5ab), 
The burst searcher software can be found at [Heimdall](https://sourceforge.net/projects/heimdall-astro/). The burst classifier software FETCH is be found [here](https://github.com/devanshkv/fetch). The pulsar package DSPSR is hosted at [here](https://sourceforge.net/projects/dspsr/) and SIGPROC can be retrieved from [here](https://github.com/SixByNine/sigproc).
    
Files:

    frb.yml:
    YAML file to reproduce the Python environment for the analysis.

    Bash scripts (.sh):
    Scripts used to generate and prepare filterbank files and archive files used in the analysis.

Folders:

    ./databases:
    The databases (.csv files) are stored in ./databases/ and are structured as follows:
        FRB20240619D_HyperFlash_paper_table.csv: Table file correspondig to burst properties for the HyperFlash detected bursts - Table D1
        FRB20240619D_ECLAT_paper_table.csv : Table file corresponding to burst properties for the ÉCLAT detected bursts - Table D2

        hyperflash_mkt_burst_db.csv - raw information on all HyperFlash burst components
        NRT_RM_values_and_errors.csv - ÉCLAT polarization values and it's errors in raw format
        /nrt/nrt_burst_db_v2.csv - raw information of NRT burst properties

    ./dm_analysis:
        /Appendix_F2_NRT_bright_burst
            Code to recreate Appendix Figure F1
            
        /Fig1_dm_over_time
            Code to recreate Figure 1 of the paper
            Also includes code of the linear fits on the DM and RM measurments

        /Fig2_Westerbork_baseband_bursts
            "Dynamix Matrix" - Figure 2 of the paper in which we show the DM evolving over time.

        /dm_opt_baseband
            Code and structure to run SFXC in a for-loop in order to create the coherent S/N vs. DM curves shows in Fig 2.
            The .csv files are present in the folders.
            
    ./python:
    Relevant Python scripts and Jupyter Notebook code used to analyze the data and generate the plots. Each file includes relevant comments.

        /Burst-rates.ipynb
          The code to create the cumulative energy distribution plot
          This code also explains all steps and sources data from .csv files in this repo
          
    ./flags:
    Flag files for each burst used for zapping radio frequency interference (RFI).

    ./plots:
    Folder containing all plots in .pdf and/or .png format.

    ./spc:
    Bash script used to execute the SPC algorithm on .archive files.

If you used this code or any of the tools, please cite [Ould-Boukattine+2026](https://arxiv.org/abs/2509.16374)  in your work.
