#!/bin/bash
set -e

##########################
#
# Once we have created all the ar-files (SFXC, scaled, SPC'ed), this
# script is wrapper around .python/ACF_archive.py which lets you compute
# the ACF, the widths, the S/N, the fluences for all burst components. All
# data are stored in a pandas data frame.
#
# Stockert command: ls *ds.ar | sort -t '/' -k 4 | cut -d '.' -f 1 | cut -c 8- |sed -e 's;^;#"st_mkt st ;' | sed '$!s/$/ ${sefd_st} B999-st" \\/'
##########################

sefd_wb=420
sefd_st=385

# 2025-03-26: First run of all bursts (without flagging of subband edges)

runs=( \
# exp dish scan sefd ID
#"p55025 wb 28 L 464.857 ${sefd_wb} B01-wb D" \
#"p55038 wb 23 L 464.857 ${sefd_wb} B02-wb D" \
#"st_mkt st 60526_8897634 L 464.857 ${sefd_st} B02-st D" \
#"st_mkt st 60530_8765030 L 464.857 ${sefd_st} B03-st D - wide - time4" \
#"st_mkt st 60535_9395728 L 464.857 ${sefd_st} B04-st D - bright" \
#"p55048 wb 24 L 464.857 ${sefd_wb} B05-wb D" \
#"st_mkt st 60536_8782036 L 464.857 ${sefd_st} B05-st X - faint - time4" \
#"p55049 wb 25 L 465.095 ${sefd_wb} B06-wb D" \
#"st_mkt st 60553_8931087 L 464.857 ${sefd_st} B07-st X - faint - time4" \
#"p55065 wb 29 L 465.294 ${sefd_wb} B07-wb D" \
#"p55074 wb 30 L 464.857 ${sefd_wb} B10-wb D" \
#"p55081 wb 19 L 465.266 ${sefd_wb} B11-wb D" \

#2 bit upsampled bursts
#"p55069 wb 27 L 464.857 ${sefd_wb} B08-wb D" \
#"p55069 wb 31 L 464.857 ${sefd_wb} B09-wb D" \

)

basedir=/data1/omar/
dbdir=/home/omar/git/frb20240619d-single-dish-analysis/dbs/

#Scrunch is done per burst, depending on visibilty of the burst. 
# for O8, Wb, Tr
Tscrunch=1
Fscrunch=1

#Distance to the source is not known at the time of analysis
distance=999.99  # Mpc

srcs=( "sfxc" "scale" "spc" )    # used as tag in the data frame to designate source of the info
compare_to=( " " " " "--compare_to scale" )

ACF_archive=./python/ACF_archive.py
ACF_archive_params="--newfit --distance ${distance} --show_acf_plots"

for run in "${runs[@]}";do
    IFS=" " read -r -a info <<< "${run}"
    exp="${info[0]}"
    dish="${info[1]}"
    scan="${info[2]}"
    band="${info[3]}"
    DM_dspsr="${info[4]}"
    sefd="${info[5]}"
    id="${info[6]}"

    arbase=${exp}_${dish}_no00${scan}

    #Setting DM
    DM_sfxc=0
    if [[ ${exp} == 'p55069' ]];then
        DM_sfxc=464.857
    fi
    
    db="${dbdir}/hyperflash_mkt_burst_db.csv"  # everything goes here
    echo "Setting dspsr DM to $DM_dspsr"
    
    dms=( ${DM_sfxc} ${DM_dspsr} ${DM_dspsr} )

    #The name of stockert burst will be the MJD of detection
    if [[ ${dish} == 'st' ]];then
        arbase=${exp}_${scan}
    fi

    expbasedir=${basedir}/${exp}
    sfxc_file=${expbasedir}/sfxc/${arbase}.ds.pazi
    scale_file=${expbasedir}/scale_spc/${arbase}_allIFs.vdif_pol2.fil.scale.ds.pazi
    spc_file=${expbasedir}/scale_spc/${arbase}_allIFs.vdif_pol2.fil.scale.spc.ds.pazi
    fs="${sfxc_file} ${scale_file} ${spc_file}"
    fscrunch=${Fscrunch}
    tscrunch=${Tscrunch}
    cnt=0

    for f in ${fs}; do
        if ! [ -e ${f} ];then
            echo "${f} doesn't exist. Moving on..."
            let cnt=${cnt}+1
            continue
        fi
        #db="${dbs[${cnt}]}"
        dm="${dms[${cnt}]}"
        src="${srcs[${cnt}]}"
        comp_to="${compare_to[${cnt}]}"
        if [[ ${dish} == 'st' ]];then
            dm=${DM_dspsr}
            #fscrunch=1
            #tscrunch=1
        fi
        cmd="python3 ${ACF_archive} ${ACF_archive_params} --id ${id} --fscrunch ${fscrunch} --tscrunch ${tscrunch} -d ${dm} --db ${db} --src ${src} -S ${sefd} ${comp_to} ${f}"

        echo "running ${cmd}"
        eval ${cmd}
        let cnt=${cnt}+1
    done
done
