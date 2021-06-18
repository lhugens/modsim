#!/bin/bash

for i in  {1..30};
PRIM=0.1;
SEC=0.156;
do
rm readingpressure.txt
touch readingpressure.txt
STR= $(($PRIM*${i})); 
STR= $(($SEC*$STR));
echo $STR;
done
