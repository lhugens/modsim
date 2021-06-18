#!/bin/bash

for k in {0..49};
do
a=$(($((2*$k))+1));
echo $a;
done
