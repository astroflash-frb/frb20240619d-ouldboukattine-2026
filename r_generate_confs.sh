#!/bin/bash
#####
# script ran on 17 March 2025
# The time and freq resolution from running base2fil does not matter for the SPC algorithm (-n 512 -d 1), since we will not be using the filterbank files, but are intrested in the VDIF frames 
# The 'keepVDIF' split .vdif frames are the important step.
# IN the base2fil script there are also .header (HDR) files written which digifil needs in the spc_ar file script, so we need to write out the .fil files nontetheless

# Before running the created configs, one has to unmount possibly mounted observations (fusermount -u home/user/vbs_data/exp)

# The skip is roughly two seconds before the burst according to digifil for L-band

# for bursts in same scan, cp the vdif file and HDR files with a 2 appended to it, this way you can process multiple bursts. 

# exp dish scan skip
#"r67008 o8 17 3  " \
####
runs=( \
"p55025 wb 28 495" \
"p55038 wb 23 127" \
"p55048 wb 24 16" \
"p55049 wb 25 258" \
"p55065 wb 29 464" \
"p55074 wb 30 323" \
"p55081 wb 19 816" \
)

#Lband 
lengths=5 

expcounter=-1
for run in "${runs[@]}";do
    IFS=" " read -r -a info <<< "${run}"
    exp="${info[0]}"
    dish="${info[1]}"
    scan="${info[2]}"
    skip="${info[3]}"
    echo "${exp} ${dish} ${scan} ${skip}"
    conf_file=./confs/${exp}_${dish}_Rmkt_scan${scan}.conf

    #if statement on source, if tr or wb
    if [ "${dish}" == "wb" ];then
	    source="R240619D"
    fi

    create_config.py -t ${dish} -s ${source} \
                     --keepVDIF --keepBP --nbit=-32 -n 512 \
                     -N 20 -i ./vex_files/${exp}.vex \
		             -S ${scan} -o ${conf_file}
    sed -i -e "s/skips.*/skips=( ${skip} )/" ${conf_file}
    sed -i -e "s/lengths.*/lengths=( ${lengths} )/" ${conf_file}
done

#    if [[ ${exp} == 'r67022' ]]; then
#	sed -i -e "s/freqLSB_0.*/freqLSB_0=1468.0/" ${conf_file}
#    fi
