#!/bin/bash
#
#########
# Ould-Boukattine 2024/2025
# Kirsten & Ould-Boukattine 2023

# This script assumes that the VDIF files were already split
# into individual subbands and live somewhere on disk. This can be done using the base2fil script.
# It performs the following steps:
# - using 'digifil', it will create Stokes I filterbank files for each subband
# - using 'dspsr', it will create archives files, attempting to have the pulse in the middle
# - each ar-file (one per subband) will first be scaled -- there can't be smearing left in a channel for scaling to work
# - the scaled ar-files with then be corrected using psrchive's SPC alogorithm
# - the scaled ar-files will be co-added to have the full band in one archive
# - the first scaled and then SPC'ed ar-files will be co-added to have the full band in one archive
# - a diff-archive of the scaled and the SPC'ed file will be generated
# - diagnostic plots of the scaled, SPC'ed, diff'ed archives are generated
#########

bw64=64000000
bw128=128000000
bw256=256000000
bw512=512000000

subbw8=8000000
subbw16=16000000
subbw32=32000000

pwait() {
    while [ $(jobs -p | wc -l) -ge $1 ]; do
	sleep 0.33
    done
}
njobs=16
expcounters=( 0 ) # this is zero based!

# You cant place a space after the \ at the eol
runs=( \
# exp - dish -scan - bw - subbandbw - cepoch_s -scaling factor for spc - index
#"p55025 wb 28 464.857 ${bw128} ${subbw16} 0.1 0.50 0 D 1sec" \
#"p55038 wb 23 464.857 ${bw128} ${subbw16} 0.2 0.50 1 D" \
#"p55048 wb 24 464.857 ${bw128} ${subbw16} 0.4 0.50 2 D 2sec&start_b1.5sec" \
#"p55049 wb 25 465.095 ${bw128} ${subbw16} 1.95 0.50 3 D" \
#"p55065 wb 29 465.294 ${bw128} ${subbw16} 1.35 0.50 4 D 2sec" \
#"p55074 wb 30 464.857 ${bw128} ${subbw16} 1.85 0.50 5 D" \
"p55081 wb 19 465.266 ${bw128} ${subbw16} 0.0 0.60 6 D" \
)

# 4.15e6 * 464.857 * (1207^-2 - (1207+0.015625)^-2) = 0.034 ms = 34 us
# We need to create 64 us with 15.625 KHz channels in order for SPC to work
# Remember, we can' have any smearing in the channels, otherwise SPC will not work
# This chosen frequency is also at Nyquist.

# The sweep at this DM is ~250 ms over the entire band
# We need to therefore make sure that we sample the entire burst

# What one can do is after applying the SPC algortithm is downsample if frequency, 
# that is what the pam command are for. Then we go back to the 125 KHz. 

# TODO implement different time- and freq- options for bursts, BUT needs to be the same for SFXC, apples to apples

######FREQ resolutions
#freq_resolution=3906.25
#freq_resolution=7812.5
freq_resolution=15625.0 
#freq_resolution=31250.0  
#freq_resolution=62500.0  
#freq_resolution=31250.0  

#######TIME RESOLUTIONS
time_res=0.000064
#time_res=0.000128
#time_res=0.000256 
#time_res=0.000512
#time_res=0.001024
#time_res=0.002048

#Downsample commands at the end of the script
pam_tscrunch=2 #to get it to 128 us
pam_fscrunch=32 #To get it to 500 KHz, the same as most of the SFXC created filterbanks
    
basedir=/data1/omar/
vdif_base=/scratch0/omar/
dir_in_outdir=scale_spc/

#Latest installs of dspsr and digifil on Sleipnir (March-2025)
dspsr=/home/oper/software/installs/dspsr-psrchive-latest/bin/dspsr
digifil=/home/oper/software/installs/dspsr-psrchive-latest/bin/digifil
pam=/home/oper/software/installs/dspsr-psrchive-latest/bin/pam

# we create archives on a per subband level, then scale them, spc them, psradd them.
# need to have NPOL 2 in hdr files for vdif-devel dspsr to work.
#for expcounter in "${expcounters[@]}";do
#    run="${runs[${expcounter}]}"
for run in "${runs[@]}";do
    IFS=" " read -r -a info <<< "${run}"
    exp="${info[0]}"
    dish="${info[1]}"
    scan="${info[2]}"
    dm="${info[3]}"
    bw="${info[4]}"
    subbandbw="${info[5]}"
    cepoch_s="${info[6]}"
    scale="${info[7]}"

    nif=`echo "${bw}/${subbandbw}" | bc | cut -d '.' -f1`
    
    nif_used=${nif}    
    echo "time res: $time_res and freq res: $freq_resolution"

    #length=2.097152
    length=1.048576

    #Lenght digifil filterbanks; account for the sweep while leaving guassian noise in the file
    digi_l=3.9

    #The .vdif frames on disk are 5 seconds long (L-band), we want the middle 2 seconds with the burst, so we Seek 1 second into the file.
    start_b=1.0

    echo "Setting DM to $dm"

    echo "${exp} ${dish} ${band} ${scan} ${bw} ${subbandbw} ${cepoch} ${nif_used} ${dm}"
    nbins=`echo "${length}/${time_res}" | bc -l | cut -d '.' -f1`

    #cepoch is given in seconds above, here we convert to MJD and add a leading 0. 
    cepoch=`echo "(( ${cepoch_s} / ( 24 * 3600 ) ))" | bc -l | sed -e 's/^-\./-0./' -e 's/^\./0./'`

    #Define the outdir
    outdir=${basedir}/${exp}/${dir_in_outdir}

    #Make the outdir in case did not exist.
    if ! [ -d ${outdir} ];then
	    mkdir -p ${outdir}
    fi
    
    nchans=`echo "${subbandbw}/${freq_resolution}" | bc | cut -d '.' -f1`
    leakage_factor=`echo "${nchans}*4" | bc | cut -d '.' -f1`

    for IF in `seq 1 ${nif_used}`;do
	    hdr=${exp}_${dish}_no00${scan}_IF${IF}.vdif_pol2.hdr

        #Remove .hdr 
	    ar=${outdir}/${hdr%.hdr}

        #Remove any filterbanks that were already there in the output folder (handy in case we rerun)
        rm ${ar}.fil

        #Digifil command
        cmd="${digifil} -I0.0 -b8 -D0.0 -d1 -F${nchans}:${leakage_factor} -T ${digi_l} -o ${ar}.fil ${vdif_base}/${exp}/${hdr} && "

        echo "ALERT ALERT ALERT THIS IS THE DM USED: $dm"
        cmd=${cmd}"${dspsr} -S ${start_b} -b ${nbins} -D ${dm} -c ${length} -T ${length} -cepoch ${cepoch} ${ar}.fil -O ${ar}.fil && "
        cmd=${cmd}"./spc/scale_${scale}.sh -e scale.ar ${ar}.fil.ar && "
        cmd=${cmd}"./spc/spc.sh -e spc.ar ${ar}.fil.scale.ar && "
        cmd=${cmd}"${pam} -f ${pam_fscrunch} -b ${pam_tscrunch} -e ds.ar ${ar}.fil.scale.ar && "
        cmd=${cmd}"${pam} -f ${pam_fscrunch} -b ${pam_tscrunch} -e ds.ar ${ar}.fil.scale.spc.ar "
        echo "$cmd"
        eval ${cmd} &
        pwait ${njobs}
    done
    wait

    # add the subbands into one
    #Define the files
    ar_scale=${exp}_${dish}_no00${scan}_allIFs.vdif_pol2.fil.scale.ar
    ar_spc=${exp}_${dish}_no00${scan}_allIFs.vdif_pol2.fil.scale.spc.ar
    ar_scale_ds=${exp}_${dish}_no00${scan}_allIFs.vdif_pol2.fil.scale.ds.ar
    ar_spc_ds=${exp}_${dish}_no00${scan}_allIFs.vdif_pol2.fil.scale.spc.ds.ar

    #Using psradd we will stitch the IFs together 
    psradd -R -o ${outdir}/${ar_scale} ${outdir}/${exp}_${dish}_no00${scan}_IF*.vdif_pol2.fil.scale.ar
    psradd -R -o ${outdir}/${ar_scale_ds} ${outdir}/${exp}_${dish}_no00${scan}_IF*.vdif_pol2.fil.scale.ds.ar
    psradd -R -o ${outdir}/${ar_spc} ${outdir}/${exp}_${dish}_no00${scan}_IF*.vdif_pol2.fil.scale.spc.ar
    psradd -R -o ${outdir}/${ar_spc_ds} ${outdir}/${exp}_${dish}_no00${scan}_IF*.vdif_pol2.fil.scale.spc.ds.ar

    #Making diag plots

    #Plot of the downsampled scaled_archived file
    psrplot -pfreq+ -c x:unit=s -D /CPS -jpDT -j 'F x1' -j 'B x8' ${outdir}/${ar_scale_ds}
    #psrplot -pfreq+ -c x:unit=s -D /CPS -jpDT -j 'F x4' -j 'B x8' ${outdir}/${ar_scale}
    mv pgplot.ps ${outdir}/${ar_scale_ds}.ps

    #Plot of the downsampled spc_archive file
    psrplot -pfreq+ -c x:unit=s -D /CPS -jpDT -j 'F x1' -j 'B x8' ${outdir}/${ar_spc_ds}
    mv pgplot.ps ${outdir}/${ar_spc_ds}.ps

    #psrdiff is to find the difference in profile; see the impact of the spc algorithm
    #Specifically the downsampled archive files
    psrdiff ${outdir}/${ar_scale_ds} ${outdir}/${ar_spc_ds}
    mv psrdiff.out ${outdir}/${exp}_${dish}_no00${scan}_scale_spc_diff_ds.ar
    psrplot -pD -c x:unit=s -D /CPS -jpFD ${outdir}/${exp}_${dish}_no00${scan}_scale_spc_diff_ds.ar
    mv pgplot.ps ${outdir}/${exp}_${dish}_no00${scan}_scale_spc_diff_ds.ar.ps

done
