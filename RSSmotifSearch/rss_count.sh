#!/bin/bash

RSS_FILE=$1
FASTA_FILE=${RSS_FILE/_rss.txt/.fasta}
NAME=$(basename ${RSS_FILE})

RSS1=($(grep "PASS$" ${RSS_FILE} | grep "^chr" | cut -f 1 | grep -P "1$" | sort | uniq ))
RSS2=($(grep "PASS$" ${RSS_FILE} | grep "^chr" | cut -f 1 | grep -P "2$" | sort | uniq ))
BREAKS=$(grep -P "^>.+1$" ${FASTA_FILE} | wc -l )

RSS_BOTH=()


for R1 in ${RSS1[@]}; do
   R1=${R1/_1/}
   for R2 in ${RSS2[@]}; do
       R2=${R2/_2/}
       if [ "${R1}" == "${R2}" ]; then
          RSS_BOTH+=(${R1})
       fi
   done
done


echo ${NAME} ${#RSS1[@]} ${#RSS2[@]} ${#RSS_BOTH[@]} ${BREAKS}