#!/bin/bash

if [ -f FILES ]; then
    rm FILES
fi
replicate=1
while [ $replicate -le 5 ]
do
    echo "10000_${replicate}_0.01.results" >> FILES
    replicate=$(( $replicate + 1 ))
done
python plots.py < FILES
if [ -f FILES ]; then
    rm FILES
fi
