#!/bin/bash

if [ -f FILES ]; then
    rm FILES
fi
replicate=1
while [ $replicate -le 5 ]
do
    m=15
    while [ $m -le 150 ]
    do
	echo "${m}_10000_${replicate}_0.01.txt" >> FILES
	m=$(( $m + 5 ))
    done
    python complexity.py "10000_${replicate}_0.01.results" < FILES
    replicate=$(( $replicate + 1 ))
    if [ -f FILES ]; then
	rm FILES
    fi
done
