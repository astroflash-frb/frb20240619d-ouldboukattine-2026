#!/bin/bash
set -e

###############
#
# This script takes the SFXC-generated, coherently dedispersed filterbanks
# and downsamples them to the same time and frequency resolution as what the
# digifil-generated filterbanks have. It then puts the files in the right
# directories.
#
# These filterbanks were used to determine the TOAs of the bursts.
#
# Generate St  command: ls *.fil | sort -t '/' -k 4 | sed -e 's;^;#"stocke st 01 01 ;' | sed '$!s/$/" \\/'
###############

sfxc_dir=/data1/omar/sfxc/frb240619d/
outbasedir=/data1/omar/

runs=( \
# exp  dish scan dir filterbank
#"p55025 wb L 28 p55025 p55025_28_1us_500KHz_StokesI_DM464857.cor_Wb.fil" \
#"p55038 wb L 23 p55038 p55038_23_128us_125KHz_StokesI_DM464857.cor_Wb.fil" \
#"p55048 wb L 24 p55048 p55048_24_32us_125KHz_StokesI_DM465106.cor_Wb.fil" \
#"p55049 wb L 25 p55049 p55049_25_1us_500KHz_StokesI_DM465095.cor_Wb.fil" \
#"p55065 wb L 29 p55065 p55065_29_64us_125KHz_StokesI_DM465294.cor_Wb.fil" \
#"p55074 wb L 30 p55074 p55074_30_128us_125KHz_StokesI_DM464857.cor_Wb.fil" \
#"p55081 wb L 19 p55081 p55081_19_1us_500KHz_StokesI_DM465266.cor_Wb.fil" \

#Stockert
"st_mkt st L 01 01 snip60526.8897634.fbspi.fil" \
"st_mkt st L 01 01 snip60530.8765030.fbspi.fil" \
"st_mkt st L 01 01 snip60535.9395728.fbspi.fil" \
"st_mkt st L 01 01 snip60536.8782036.fbspi.fil" \
"st_mkt st L 01 01 snip60553.8931087.fbspi.fil" \
)

for run in "${runs[@]}";do
    IFS=" " read -r -a info <<< "${run}"
    exp="${info[0]}"
    dish="${info[1]}"
    band="${info[2]}"
    scan="${info[3]}"
    dir="${info[4]}"
    sfxc_fil="${info[5]}"
    fulldir=${sfxc_dir}/${dir}_${scan}
    if [[ ${dish} == 'st' ]];then
        fulldir=/data1/omar/frb240619d/stockert
    fi
    target=${fulldir}/${sfxc_fil}
    if ! [[ -e ${target} ]];then
        echo "cannot find ${target}"
    fi

    #Settings of output filterbanks
    if [[ "$band" = "L" ]]; then
    ##Lband
    #This setting also works for Stockert, since it won't actually downsample
        tres_out=128         # in microseconds
        fres_out=500000     # in Hz
    elif [[ "$band" = "P" ]]; then
        tres_out=512       # in microseconds
        fres_out=31250     # in Hz
    fi

    #The name of stockert burst will be the MJD of detection
    if [[ ${dish} == 'st' ]];then
        #The name of stockert burst will be the MJD of detection
        #snip59867.5457019.fbspi.fil
        #snip59867_5457019.fbspi.fil -- replace the first dot
        sfxc_fil_d=${sfxc_fil/./_}

        #echo snip59949_8819706.fbspi.fil | cut -d '.' -f 1 | cut -c 5-
        #return 59949_8819706
        mjd_dir=`echo ${sfxc_fil_d} | cut -d '.' -f 1 | cut -c 5-`
        linkname=${exp}_${mjd_dir}.fil
    else
        linkname=${exp}_${dish}_no00${scan}.fil
    fi

    linkdir=${outbasedir}/${exp}/sfxc/
    link=${linkdir}/${linkname}

    if [ -e ${link} ];then
        rm ${link}
    fi
    if ! [ -d ${linkdir} ];then
        mkdir -p ${linkdir}
    fi
    cmd="ln -s ${target} ${link}"
    echo "running ${cmd}"
    eval ${cmd}

    ##
    #pars=`header ${link} -tsamp -nifs -foff -nbits`
    #echo ${pars}
    #IFS=" " read -r -a params <<< "${pars}"
    tres_in=`header ${link} -tsamp`   # in microseconds
    npol=`header ${link} -nifs`    # 4 for full pol, we want Stokes I
    fres_in=`header ${link} -foff`
    nbits=`header ${link} -nbits`
    if [[ ${nbits} == '32' ]];then
        nbits='-32'
    fi
    ## bash cannot handle floats, so convert fres_in to an integer in Hz
    fres_in=`echo "${fres_in}*1000000" | bc | cut -d '.' -f1`
    if [[ ${fres_in} -lt 0 ]]; then
        fres_in=`echo ${fres_in}*-1 | bc`
    fi
    ## get the scrunching factors
    tscale=`echo "${tres_out}/${tres_in}" | bc`
    fscale=`echo "${fres_out}/${fres_in}" | bc`
    echo "tres_in=${tres_in}, fres_in=${fres_in}"
    echo "tscale=${tscale}, fscale=${fscale}"
    cmd="digifil ${link} -o ${link%fil}ds.fil -I0 -b ${nbits} "
    if [[ ${tscale} -gt 1 ]]; then
        cmd=${cmd}"-t ${tscale} "
    fi
    if [[ ${fscale} -gt 1 ]]; then
        cmd=${cmd}"-f ${fscale} "
    fi
    if [[ ${npol} -gt 1 ]];then
        cmd=${cmd}"-d1 "
    fi
    echo "running ${cmd}"
    ##exit 1
    eval ${cmd}
done
