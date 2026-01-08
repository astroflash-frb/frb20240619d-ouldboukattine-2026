#!/bin/bash

# Ould-Boukattine (2025)
# Kirsten & Ould-Boukattine (2023) 
# Script takes the cut out filterbank files as provided by Wolfgang and converts them to archive files with the burst in the middle

# The time sampling of Stockert: 218.45333 us, we create the archive file to be 2048 bins
# This will yield an archive file with lenght: 0.447391744 seconds

# The script also checks if the final time resoltuon is a factor of 2, if this is not the case (7.999us-bins) for example
# Then this can be solved by changing the start_time of the archive file (by 0,01), it's behavior I suspect has to do with rounding of start time

# dont forget the extra space at the end of the line, before the \, otherwise it won't run
# middle = location of the burst in the sfxc filterbank, assumed that the burst is at the 0.5 seconds mark
# cepoch = shift to left in seconds

# logs
# 

runs=( \
# middle cepoch_s filterbank
# L-band
#"1.0 0.42 snip60526.8897634.fbspi.fil" \
#"1.0 0.22 snip60530.8765030.fbspi.fil" \
#"1.0 0.08 snip60535.9395728.fbspi.fil" \
#"1.0 0.42 snip60536.8782036.fbspi.fil" \
"1.0 0.1 snip60553.8931087.fbspi.fil" \
)

#Create directory for the archive files
basedir=/data1/omar/sfxc/frb240619d
archive_dir=/fil_to_archives_stockert

fil_dir=/data1/omar/frb240619d/stockert/

outdir_archive=${basedir}/${archive_dir}
if ! [ -d ${outdir_archive} ];then
    mkdir -p ${outdir_archive}
fi

#Define constants
dm=464.857
#218.453333333333
time_res=0.000218453 #time_resolution of the stockert recording
lenght_archive=0.447391744 #fixed, 2048 bins, see above

#calculate the number of bins based on the lenght and time_res
nbins=`echo "${lenght_archive}/${time_res}" | bc -l | cut -d '.' -f1`

for run in "${runs[@]}";do

    #Read the run line using the read function
    read -r -a info <<< "${run}"
    middle_b="${info[0]}"
    cepoch_s="${info[1]}"
    filterbank="${info[2]}"

    #Start time of the archive file
    start_b=`echo "${middle_b} - ${lenght_archive}/2" | bc -l`

    #From documentation of CEPOCH
    #To specify a different epoch (e.g. the middle of the observation) use the EPOCH key with value specified as MJD.
    #So we convert from seconds to MJD to define the shift to left.
    #In order for CEPOCH to work you need a 0. infront of the number, so using sed we add a 0, see link below
    cepoch_mjd=`echo "(( ${cepoch_s} / ( 24 * 3600 ) ))" | bc -l | sed -e 's/^-\./-0./' -e 's/^\./0./'`

    echo "${exp} ${dish} ${scan} ${direc} ${cepoch_mjd} ${filterbank}"

    filterbank_in_dir=${fil_dir}/${filterbank}
    archive_out_nm=${outdir_archive}/${filterbank/./_}
    echo $archive_out_nm

    #runnig dspsr to create an archive file
    #need to set machine as fake, otherwise we get an error regarding the processing machine pulsar2000 having to long a lenght
    dspsr_cmd="dspsr -set machine=fake -D ${dm} -c ${lenght_archive} -T ${lenght_archive} -S ${start_b} -b ${nbins} -cepoch ${cepoch_mjd} -O ${archive_out_nm} ${filterbank_in_dir}"
    echo $dspsr_cmd
    ${dspsr_cmd}

    #Name of the archive file with extension
    archive_out_nm_ar=${archive_out_nm}.ar

    #In case the RFI is so bad and you need to zap channels in order to see them in the diagnostic plot
    #First use the Pazi tool to figure our the channels, add them here below between "" and the script will zap them and the diagnostic plot will be updated accordinly.
    #paz -z "94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138" -e pazi ${archive_out_nm_ar}
    #archive_out_nm_ar=${archive_out_nm}.pazi

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
        continue
    fi

    #adding the info regarding the DM, machine and DMC to the archive file
    #Add the DM value 
    #dspsr_cmd_dm="psredit -m -c dm=${dm} ${archive_out_nm_ar}"
    #echo $dspsr_cmd_dm
    #${dspsr_cmd_dm}

    #Add the fact that the burst is already straight, dmc=1
    #dspsr_cmd_dmc="psredit -m -c dmc=1 ${archive_out_nm_ar}"
    #echo $dspsr_cmd_dmc
    #${dspsr_cmd_dmc}
     
    #Downsample factor for diagnostic plots
    time_down=2
    freq_down=1
    
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
    #echo "${new_name_tf}"
    mv pgplot.ps ${new_name_tf}.ar.ps

done

#Regarding the sed expression
#https://stackoverflow.com/questions/8402181/how-do-i-get-bc1-to-print-the-leading-zero

#making into a ps plot
#psrplot -pfreq+ -c x:unit=s -D /CPS -jpDT -j 'F x4' -j 'B x8' ${outdir}/${ar_scale}
#mv pgplot.ps ${outdir}/${ar_scale}.ps

#Selecting all filterbanks for in this file
#ls *.fil | sort -t '/' -k 4 | sed -e 's;^;#"0.5 0.0 ;' | sed '$!s/$/" \\/'

#pazi command:
#snip59870.7739324 -- paz -z "143 144 145 146 147 141 142 148 149 140 150 151 101 102 103 139" 
