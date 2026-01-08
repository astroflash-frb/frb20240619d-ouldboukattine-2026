#!/bin/bash

#####################
#
# This script takes the flag files 
# into the git repo and applies them to the scaled, SPC'ed, and SFXC archives with 'paz'.
#
# For stockert we automatically make a zap file with 0 and 168; if that is not enough I used paz to zap the channels manually and add them to the flag file
#
# Once all the flags are applied, we can start to compute fluences, S/N, widths, TOAs...
#
# Command for stockert list
# ls *ds.ar | sort -t '/' -k 4 | cut -d '.' -f 1 | cut -c 8- |sed -e 's;^;#"st_mkt st ;' | sed '$!s/$/" \\/'
# First prep, then flag
####################


pwait() {
    while [ $(jobs -p | wc -l) -ge $1 ]; do
	sleep 0.33
    done
}
njobs=16

outbasedir=/data1/omar/
flagrundir=/home/omar/git/frb20240619d-single-dish-analysis/flags/

runs=( \
# exp dish scan flag
#"p55025 wb 28" \
#"p55038 wb 23"  \
#"p55048 wb 24" \
#"p55049 wb 25" \
#"p55065 wb 29" \
#"p55074 wb 30" \
#"p55081 wb 19" \
"p55069 wb 27"
"p55069 wb 31"
#"st_mkt st 60526_8897634" \
#"st_mkt st 60530_8765030" \
#"st_mkt st 60535_9395728" \
#"st_mkt st 60536_8782036" \
#"st_mkt st 60553_8931087" \
)

for run in "${runs[@]}";do
    IFS=" " read -r -a info <<< "${run}"
    exp="${info[0]}"
    dish="${info[1]}"
    scan="${info[2]}"
    band="${info[3]}"
    arbase=${exp}_${dish}_no00${scan}
    flagfile=${arbase}.sh
    echo ${flagfile}

    #The name of stockert burst will be the MJD of detection
    if [[ ${dish} == 'st' ]];then
        arbase=${exp}_${scan}
    fi

    expbasedir=${outbasedir}/${exp}
    if [ "$dish" = "st" ] && ! [ -e ${flagrundir}/${flagfile} ]; then
        # this is for Stockert only.
        echo -e '#!/bin/bash\npaz -e pazi -z "0 167" $1' > ${flagrundir}/${flagfile}
        echo "created ${expbasedir}/${flagfile} containing"
        cat ${flagrundir}/${flagfile}
    fi

    if [ "$dish" = "wb" ] && ! [ -e ${flagrundir}/${flagfile} ]; then
        # this is for Stockert only.
        echo -e '#!/bin/bash\npaz -e pazi -z "0" $1' > ${flagrundir}/${flagfile}
        #echo -e '#!/bin/bash\npaz -e pazi -z "0,1,254,255,256,257,510,511,512,513,766,767,768,769,1022,1023,1024,1025,1278,1279,1280,1281,1534,1535,1536,1537,1790,1791" $1' > ${flagrundir}/${flagfile}
        echo "created ${expbasedir}/${flagfile} containing"
        cat ${flagrundir}/${flagfile}
    fi

    chmod u+x ${flagrundir}/${flagfile}
    sfxc_file=${expbasedir}/sfxc/${arbase}.ds.ar
    scale_file=${expbasedir}/scale_spc/${arbase}_allIFs.vdif_pol2.fil.scale.ds.ar
    spc_file=${expbasedir}/scale_spc/${arbase}_allIFs.vdif_pol2.fil.scale.spc.ds.ar
    #fs="${scale_file} ${spc_file}"
    fs="${sfxc_file} ${scale_file} ${spc_file}"
    for f in ${fs}; do
        if ! [ -e ${f} ];then
            echo "${f} doesn't exist. Moving on..."
            continue
        fi
        cmd="${flagrundir}/${flagfile} ${f}"
        echo "running ${cmd}"
        eval ${cmd}
    done
done

#pazi command:
#When creating the archives files I had to zap channels otherwise the bursts were not visible. First I let this script create the .sh files, then added these channels to that file to zap the bad channels
# p11739_wb_no0018.ds.pazi paz -z "521 522 523 524 525 526 527 528 529 530 531 532 533 534 535 536 537 538 539 540 541 542 543 544 545 546 547 548 549 550 551 552 553 554 555 556 557 558 559 560 561 562 563 564 565 566 567 568 569 570 571 572 573 574 575 576 577 578 579 580 581 582 583 584 585 586 587 588 589 590 591 894 895 896 897" -e pazi p11739_wb_no0018.ds.pazi
#paz -z "886 887 888 889 497 498 415 416" -e pazi p11748_wb_no0020_allIFs.vdif_pol2.fil.scale.ds.ar
#513 514 515 516 517 518 519 520 521 522 41 42 43 44" -e pazi p11739_wb_no0018_allIFs.vdif_pol2.fil.scale.ds.pazi

######
#snip59870.7739324 -- paz -z "143 144 145 146 147 141 142 148 149 140 150 151 101 102 103 139" 
#snip59874.0936399 -- paz -z "92 93 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138" 
#snip59876.6610997 -- paz -z "134 135 136 137 138 94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133"
#snip59878.3813848 -- paz -z "135 136 137 138 132 133 134 108 109 110 111 94 95 96 97 126 115 116 117 100 124 125 123 105 106 130 131 98 99"
#snip59883.7188223 -- paz -z "135 136 137 138 131 132 133 134 130 104 105 106 107 94 95 96 97 98 99 100 101 108 109 110 111 112 113 114 115 116 117 102 103 127 128 129 122 123 125 126 118 119 120 121"
#snip59889.7032667 -- paz -z "135 136 137 138 129 130 131 125 126 127 128 114 115 116 117 105 106 107 94 95 96 102 103 118 119 120 132 133 112 113 98 99 101 108 122 123"
#snip59897.4974988 -- paz -z "97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 52 53 54 55 56 69 70 71 72 73 74 75 135 136 137 133 134 123 124 125 126 127 128 122 131 132"
#snip59899.8813157 -- paz -z "97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 95 96"
#snip59903_2531485 -- paz -z "94 95 96 97 98 99 100 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138"
