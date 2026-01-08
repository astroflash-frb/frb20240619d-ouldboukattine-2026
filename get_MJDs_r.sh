#!/bin/bash
set -e

##############
#
# This is a wrapper around ./python/get_MJD.py which
# combines the info from the ar-files and the pandas data
# frames to fit the bursts in the SFXC/Stockert-filterbanks with
# Guassians to determine the TOAs. For SFXC-filterbanks these
# are geocentric TOAs, while for Stockert they are topocentric at Stockert.
#
# First use prep_sfxc_filterbank to prepare the symlinks for the filterbank files
# The archive files are needed to copy the masks from the archive onto the filterbank files
#############

## References frequencies
#Middle of top subband for EVN dishes and SFXC filterbanks 
#Pband (setup P - 356-4 MHz) - SFXC 
ref_freq_wb_p=352

#Wb L-band (setup L18 - 1335-8 MHz) - SFXC 
ref_freq_wb_1g=1327
 
#St L-band (script will read the header file of the filterbank and take the Highchannel as ref freq)
ref_freq_st_1g='none'

fscrunch=1
tscrunch=1

# exp dish scan sefd ID ref_freq DM
runs=( \
#EVN dishes L-band
#"p55025 wb 28 coe B01-wb ${ref_freq_wb_1g} 0.0 D" \
#"p55038 wb 23 coe B02-wb ${ref_freq_wb_1g} 0.0 D" \
#"p55048 wb 24 coe B05-wb ${ref_freq_wb_1g} 0.0 D" \
#"p55049 wb 25 coe B06-wb ${ref_freq_wb_1g} 0.0 D" \
#"p55065 wb 29 coe B07-wb ${ref_freq_wb_1g} 0.0 D" \
#"p55074 wb 30 coe B10-wb ${ref_freq_wb_1g} 0.0 D" \
#"p55081 wb 19 coe B11-wb ${ref_freq_wb_1g} 0.0 D" \

#"st_mkt st 60526_8897634 st B02-st ${ref_freq_st_1g} 464.857 D" \
#"st_mkt st 60530_8765030 st B03-st ${ref_freq_st_1g} 464.857 D" \
#"st_mkt st 60535_9395728 st B04-st ${ref_freq_st_1g} 464.857 D" \
#"st_mkt st 60536_8782036 st B05-st ${ref_freq_st_1g} 464.857 D" \
#"st_mkt st 60553_8931087 st B07-st ${ref_freq_st_1g} 464.857 D" \

#2 bit upsampled bursts
#Doing this manually probably
)

basedir=/data1/omar/
dbdir=/home/omar/git/frb20240619d-single-dish-analysis/dbs/

get_MJD=./python/get_MJD.py

for run in "${runs[@]}";do
    IFS=" " read -r -a info <<< "${run}"
    exp="${info[0]}"
    dish="${info[1]}"
    scan="${info[2]}"
    loc="${info[3]}"
    id="${info[4]}"
    refFreq="${info[5]}"
    dm="${info[6]}"
    src='sfxc'
    src_dir='sfxc'
    freqref="--ref_freq ${refFreq}"
    if [[ ${refFreq} == 'none' ]];then
        freqref=''
    fi

    #The name of stockert burst will be the MJD of detection
    if [[ ${dish} == 'st' ]];then
        arbase=${exp}_${scan}
    else
        arbase=${exp}_${dish}_no00${scan}
    fi

    expbasedir=${basedir}/${exp}
    arfile=${expbasedir}/${src_dir}/${arbase}.ds.pazi
    filfile=${expbasedir}/${src_dir}/${arbase}.ds.fil

    #rule for burst_time 2.0 and 0.5 seconds for P- and L-band
    if [[ "$refFreq" = "$ref_freq_wb_p" ]]; then
        burst_time_s=2.0 
        db="${dbdir}/hyperflash_mkt_burst_db.csv"  # everything goes here for Pband
    elif [[ "$refFreq" = "$ref_freq_st_1g" ]]; then
        burst_time_s=1.5 
        db="${dbdir}/hyperflash_mkt_burst_db.csv"  # everything goes here for Lband
    else 
        burst_time_s=0.5
        db="${dbdir}/hyperflash_mkt_burst_db.csv"  # everything goes here for Lband
    fi

    cmd="python3 ${get_MJD} -d ${dm} --db ${db} --src ${src} --arfile ${arfile} --fscrunch ${fscrunch} --tscrunch ${tscrunch} --burst_time ${burst_time_s} --ID ${id} ${freqref} --location ${loc} ${filfile} --show_dynspec"

    echo "running ${cmd}"
    eval ${cmd}
done
