#!/bin/bash

###################
# These dataproducts are used in the determination of the burst properties. 
# Note: the other script (prep_filterbanks) are used for the ToA
# !! First prep, then flag
#
# Before runnning this script you need to:
# - Create the archive files using r117_generate-archive.sh
# - Optional: Create the flag files and store them in flags -> but you can base the flag on the spc/scale version, since we will apply the same flag to all 3 files
#
# This script takes the archive-files as generated from the
# coherently dedispersed SFXC-filterbanks and performs the
# following steps:
# - downsample in time to 64 us
# - downsample in freq to 62.5 kHz
# - generate Stokes I
# This time frequency resolution is chosen, since to be identical to the SPC/SCALE versions
#
# Command used to generate Westerbork list (made with chatgpt)
# ls *.ar | sort -t '/' -k 4 | awk -F'[_]' '{print "#\"" $1, "wb", $2, "L", $1"_"$2, $0 "\" \\"}'
#
# Command used to generate stockert list
# ls *.ar | sort -t '/' -k 4 | sed -e 's;^;#"stocke st 01 L 01 ;' | sed '$!s/$/" \\/'
###################

sfxc_dir=/data1/omar/sfxc/frb240619d/fil_to_archives
sfxc_dir_st=/data1/omar/sfxc/frb240619d/fil_to_archives_stockert
outbasedir=/data1/omar/

runs=( \
# exp dish scan sfxc
# L-band
#"p55025 wb 28 L p55025_28 p55025_28_1us_500KHz_StokesI_DM464857_cor_Wb.fil.ar" \
#"p55038 wb 23 L p55038_23 p55038_23_128us_125KHz_StokesI_DM464857_cor_Wb.fil.ar" \
#"p55048 wb 24 L p55048_24 p55048_24_128us_125KHz_StokesI_DM464857_cor_Wb.fil.ar" \
#"p55049 wb 25 L p55049_25 p55049_25_1us_500KHz_StokesI_DM465095_cor_Wb.fil.ar" \
#"p55065 wb 29 L p55065_29 p55065_29_64us_125KHz_StokesI_DM465294_cor_Wb.fil.ar" \
#"p55074 wb 30 L p55074_30 p55074_30_128us_125KHz_StokesI_DM464857_cor_Wb.fil.ar" \
#"p55081 wb 19 L p55081_19 p55081_19_1us_500KHz_StokesI_DM465266_cor_Wb.fil.ar" \
#7
#Stockert
#"st_mkt st 01 L 01 snip60526_8897634.fbspi.fil.ar" \
#"st_mkt st 01 L 01 snip60530_8765030.fbspi.fil.ar" \
#"st_mkt st 01 L 01 snip60535_9395728.fbspi.fil.ar" \
#"st_mkt st 01 L 01 snip60536_8782036.fbspi.fil.ar" \
#"st_mkt st 01 L 01 snip60553_8931087.fbspi.fil.ar" \
#2bit
"p55069 wb 27 L p55069_27 p55069_wb_no0027_IFall_vdif_pol2_8bit_10sec_snippet_fil.ar" \
"p55069 wb 31 L p55069_31 p55069_wb_no0031_IFall_vdif_pol2_8bit_10sec_snippet_fil.ar" \
)

for run in "${runs[@]}";do

    IFS=" " read -r -a info <<< "${run}"
    exp="${info[0]}"
    dish="${info[1]}"
    scan="${info[2]}"
    band="${info[3]}"
    echo "scan: $scan"
    dir_folder="${info[4]}"
    sfxc_ar="${info[5]}"

    target=${sfxc_dir}/${sfxc_ar}    
    linkdir=${outbasedir}/${exp}/sfxc/
    linkname=${exp}_${dish}_no00${scan}.ar

    echo "${target}"
    if [[ "$band" = "L" ]]; then
        tres=64        # in microseconds
        fres=500000
        #fres=62500     # in Hz
    elif [[ "$band" = "P" ]]; then
        tres=512       # in microseconds
        fres=31250     # in Hz
    fi

    #The name of stockert burst will be the MJD of detection
    if [[ ${dish} == 'st' ]];then
        #echo snip59949_8819706.fbspi.fil.ar | cut -d '.' -f 1 | cut -c 5-
        #return 59949_8819706
        #The name of stockert burst will be the MJD of detection
        mjd_dir=`echo ${sfxc_ar} | cut -d '.' -f 1 | cut -c 5-`
        linkname=${exp}_${mjd_dir}.ar
        
        #Change the location of the archives for stockert bursts
        target=${sfxc_dir_st}/${sfxc_ar}
    fi

    link=${linkdir}/${linkname}
    echo "${link}"
    if [ -e ${link} ];then
        rm ${link}
    fi
    if ! [ -d ${linkdir} ];then
	mkdir -p ${linkdir}
    fi
    ln -s ${target} ${link}
    # use $(( $var )) to remove white spaces around integers
    pars=`psredit -q -Q -c nbin,length,npol,nchan,bw ${link}`
    IFS=" " read -r -a params <<< "${pars}"
    nbin="${params[0]}"
    length="${params[1]}"  # comes in seconds
    npol="${params[2]}"
    nchan="${params[3]}"
    bw="${params[4]}" # in MHz, -ve
    # bash cannot handle floats, so convert bw to an integer in Hz
    bw=`echo "${bw}*1000000" | bc | cut -d '.' -f1`
    if [[ ${bw} -lt 0 ]]; then
        bw=`echo ${bw}*-1 | bc`
    fi
    # convert length to microseconds
    length=`echo "${length}*1000000" | bc`
    time_res=`echo "${length}/${nbin}" | bc`
    freq_res=`echo "${bw}/${nchan}" | bc`
    tscale=`echo "${tres}/${time_res}" | bc`
    fscale=`echo "${fres}/${freq_res}" | bc`
    echo "length=${length}, nbin=${nbin}, npol=${npol}"
    echo "nchan=${nchan}, bw=${bw}"
    echo "time_res=${time_res}, freq_res=${freq_res}"
    echo "tscale=${tscale}, fscale=${fscale}"
    cmd='pam -e ds.ar '
    if [[ ${tscale} -gt 1 ]]; then
        cmd=${cmd}"-b ${tscale} "
    fi
    if [[ ${fscale} -gt 1 ]]; then
        cmd=${cmd}"-f ${fscale} "
    fi
    if [[ ${npol} -gt 1 ]];then
        cmd=${cmd}"-p "
    fi
    cmd=${cmd}${link}
    echo "running ${cmd}"
    #exit 1
    eval ${cmd}
done
