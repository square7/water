#!/bin/bash

cp params.cpy params.txt
for i in `seq 1000 100 18000`
do
echo $i
python gen.py $i
./a.out > params.tmp
mv params.tmp params.txt
cat params.txt
done