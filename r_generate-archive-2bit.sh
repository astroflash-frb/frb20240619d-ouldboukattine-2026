#!/bin/bash

# Ould-Boukattine (2025)
# Kirsten & Ould-Boukattine (2023) 

# Script takes 2bit upsampled filterbanks and converts them to archive files
# This script is 95% based on the r_generate-archive script, which uses SFXC created filterbank as base

# The script also checks if the final time resoltuon is a factor of 2, if this is not the case (7.999us-bins) for example
# Then this can be solved by changing the start_time of the archive file (by 0,01), it's behavior I suspect has to do with rounding of start time

# dont forget the extra space at the end of the line, before the \, otherwise it won't run
# middle = location of the burst in the sfxc filterbank, used prepdata to see the location (plots in evernote)
# cepoch = shift to left in seconds

runs=( \
# exp   dish scan band DM direc middle cepoch_s filterbank
# 2bit upscaled L-band -- MAKE SURE TO CHANGE THE BASEDIR
#"p55069 wb 27 L 464.857 p55069_27 5.0 0.09 p55069_wb_no0027_IFall_vdif_pol2_8bit_10sec_snippet.fil D" \
"p55069 wb 31 L 464.857 p55069_31 5.0 0.03 p55069_wb_no0031_IFall_vdif_pol2_8bit_10sec_snippet.fil D" \
)

#Create directory for the archive files
basedir=/data1/omar/sfxc/frb240619d
archive_dir=/fil_to_archives

fil_dir=/data1/omar/frb240619d/wb-2bit


outdir_archive=${basedir}/${archive_dir}

if ! [ -d ${outdir_archive} ];then
    mkdir -p ${outdir_archive}
fi

for run in "${runs[@]}";do

    #Read the run line using the read function
    read -r -a info <<< "${run}"
    exp="${info[0]}"
    dish="${info[1]}"
    scan="${info[2]}"
    band="${info[3]}"
    dm="${info[4]}"
    direc="${info[5]}"
    middle_b="${info[6]}"
    cepoch_s="${info[7]}"
    filterbank="${info[8]}"

    lenght_archive=0.524288 #2**16*8e-6
        
    #For the 2bit bursts, the time resolution is 256 us.
    time_res=0.000256
        
    #naming convention:
    #filterbank_in_dir=${basedir}/${direc}/${filterbank}
    filterbank_in_dir=${fil_dir}/${filterbank}
    archive_out_nm=${outdir_archive}/${filterbank/./_}
    echo $archive_out_nm

    #calculate the number of bins based on the lenght and time_res
    nbins=`echo "${lenght_archive}/${time_res}" | bc -l | cut -d '.' -f1`

    #Start time of the archive file
    start_b=`echo "${middle_b} - ${lenght_archive}/2" | bc -l`

    #From documentation of CEPOCH
    #To specify a different epoch (e.g. the middle of the observation) use the EPOCH key with value specified as MJD.
    #So we convert from seconds to MJD to define the shift to left.
    #In order for CEPOCH to work you need a 0. infront of the number, so using sed we add a 0, see link below
    cepoch_mjd=`echo "(( ${cepoch_s} / ( 24 * 3600 ) ))" | bc -l | sed -e 's/^-\./-0./' -e 's/^\./0./'`

    echo "${exp} ${dish} ${scan} ${band} ${direc} ${cepoch_mjd} ${filterbank} ${nbins}"

    #runnig dspsr to create an archive file
    dspsr_cmd="dspsr -c ${lenght_archive} -T ${lenght_archive} -S ${start_b} -b ${nbins} -D ${dm} -cepoch ${cepoch_mjd} -O ${archive_out_nm} ${filterbank_in_dir}"
    echo $dspsr_cmd
    ${dspsr_cmd}

    #Name of the archive file with extension
    archive_out_nm_ar=${archive_out_nm}.ar

    #A check to see if the the amount of bins and output time is correct
    pars=`psredit -q -Q -c nbin,length,npol,nchan,bw ${archive_out_nm_ar}`
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
    echo "length=${length}, nbin=${nbin}, npol=${npol}"
    echo "nchan=${nchan}, bw=${bw}"
    echo "time_res=${time_res}, freq_res=${freq_res}"

    if [[ $((${time_res} % 2)) != 0 ]]; then
        echo "ALERT ALERT ALERT LOOK AT THE START TIME, current time res is: $time_res"
    fi

    ##adding the info regarding the DM, machine and DMC to the archive file
    ##Add the DM value 
    #dspsr_cmd_dm="psredit -m -c dm=${dm} ${archive_out_nm_ar}"
    ##echo $dspsr_cmd_dm
    #${dspsr_cmd_dm}

    ##Add the fact that the burst is already straight, dmc=1
    #dspsr_cmd_dmc="psredit -m -c dmc=1 ${archive_out_nm_ar}"
    ##echo $dspsr_cmd_dmc
    #${dspsr_cmd_dmc}
     
    #Downsample factor for diagnostic plots
    time_down=8
    freq_down=4
    
    freq_res_khz=`echo "${freq_res}/1000" | bc` #convert to khz

    #Generate the new time- and freq resoltuon
    plt_time_res=`echo "$time_res*$time_down" | bc` #us
    plt_freq_res=`echo "$freq_res_khz*$freq_down" | bc` #khz
    echo "Making a diagnostic plot with time_res=${plt_time_res} us and freq_res=${plt_freq_res} hz"

    psrplot -pfreq+ -c x:unit=s -D /CPS -jpDT -j 'F x'"${freq_down}" -j 'B x'"${time_down}" ${archive_out_nm_ar}

    #${parameter//pattern/string}
    #Replace the old time and freq-res with the new in the archive string
    #Only search for time_res+us, otherwise you will also replace numbers in the scan
    new_name_t=`echo "${archive_out_nm_ar/${time_res}"us"/${plt_time_res}"us"}"`
    new_name_tf=`echo "${new_name_t/${freq_res_khz}"KHz"/${plt_freq_res}"KHz"}"`
    echo "${new_name_tf}"
    mv pgplot.ps ${new_name_tf}.ar.ps

done

#Regarding the sed expression
#https://stackoverflow.com/questions/8402181/how-do-i-get-bc1-to-print-the-leading-zero

#making into a ps plot
#psrplot -pfreq+ -c x:unit=s -D /CPS -jpDT -j 'F x4' -j 'B x8' ${outdir}/${ar_scale}
#mv pgplot.ps ${outdir}/${ar_scale}.ps



