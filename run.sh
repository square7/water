#!/bin/bash

#cp params.cpy params.txt
touch ene.txt
rm ene.txt
touch ene.txt
for i in `seq 20000 5 35000`
do
echo $i
python gen.py $i
#./a.out
./a.out >> ene.txt
#./a.out > params.tmp
#mv params.tmp params.txt
#cat params.txt
done