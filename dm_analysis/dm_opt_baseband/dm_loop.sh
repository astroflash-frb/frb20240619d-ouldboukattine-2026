#!/bin/bash

#p55038_23
#arr=( $(seq 464.5 0.01 466.5) )

#p55038_23 -- notch -- expected at 464.95 
arr=( $(seq 465.6 0.005 465.9) )

#p55081_19
#arr=( $(seq 465.1 0.001 465.4) )

#p55065_29
#arr=( $(seq 465.0 0.002 465.65) )
#arr=( $(seq 464.76 1.0 464.96) )

#p55048_24
#arr=( $(seq 464.95 0.001 465.25) )

for i in "${arr[@]}";do

    #define some variables
    obscode="p55038"
    obscode_scan="${obscode}_23"
    time_res="256us"
    freq_res="500KHz"
    timestamp="0.7228"
    
    #Define some files, could be optimized somehow
    #Change these variabeles
    gitfolder="/home/omar/git/frb20240619d-single-dish-analysis/dm_analysis/dm_opt_baseband/${obscode_scan}/"
    CTRLFILE="${gitfolder}${obscode_scan}_${time_res}.ctrl"
    VIXFILE="${gitfolder}${obscode}.vix"
    
    #replace the comma with a dot, store in variable
    #echo $DMVALDOT
    DMVALDOT=${i/,/.}
    
    #find (before the c\) and replace the DM in the polyco file
    folderpolyco="/data1/omar/sfxc/frb240619d/"
    POLYCO="${folderpolyco}frb240619d_loop.polyco"
    sed -i '/frb/ c\frb20240619d   03-May-21    191325.00      59337.80   '$DMVALDOT' 0.000  0.000' "$POLYCO"

    DMVALSPACE=${i/./}
    #echo $DMVALSPACE

    #find and replace the cor name filename of the core file
    CORNAME="${obscode_scan}_${time_res}_${freq_res}_StokesI_DM$DMVALSPACE.cor"
    #echo $CORNAME

    #CORNAMECTRL="file:///dev/shm//$CORNAME"
    CORNAMECTRL="file:///data1/omar/sfxc/frb240619d/${obscode_scan}/dm_opt/$CORNAME"
    #echo $CORNAMECTRL
    sed -i '/output/ c\    "output_file": "'$CORNAMECTRL'",' "$CTRLFILE"

    #running the SFXC command
    SFXC="mpirun -np 38 --use-hwthread-cpus sfxc $CTRLFILE $VIXFILE"
    echo $SFXC
    #Excecute command
    ${SFXC}

    #Add _wb to the corname, to fit the output of sfxc. need to do a trick to get _wb appended. So use + and then replace with _
    CORNAMEOUT="${CORNAME}_Wb"
    #echo $CORNAMEOUT

    #The output filterbank name with formatting
    FILNAME="${CORNAMEOUT/./_}.fil"
    #echo $FILNAME

    #Running cor2filterbank
    folderoutput="/data1/omar/sfxc/frb240619d/${obscode_scan}/dm_opt/"
    #folderoutput="/dev/shm/"
    COR2="cor2filterbank.py -s Wb -p I $VIXFILE $folderoutput$CORNAMEOUT $folderoutput$FILNAME"
    echo $COR2
    ${COR2}

    #Calling the python script
    PYTCOM="python standard_profile_loop_min.py -f $folderoutput$FILNAME -dm $DMVALDOT -t $timestamp"
    echo $PYTCOM
    ${PYTCOM}

    #Removing the .cor and the .fil file
    RMCOM="rm $folderoutput$FILNAME $folderoutput$CORNAMEOUT"
    echo $RMCOM
    ${RMCOM}
    
done

RMCOMCHEX="rm *dynamic_channel* chex*"
echo $RMCOMCHEX
${RMCOMCHEX}
