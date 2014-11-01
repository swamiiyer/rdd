#!/bin/bash

m=15
while [ $m -le 150 ]
do
    python gen_data.py -m $m -r 5 -n 10000 -t 0.01
    m=$(( $m + 5 ))
done
