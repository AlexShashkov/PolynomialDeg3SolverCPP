#!/bin/sh

A1=$[`wc -l < $1`-1];
B1=`awk '{s+=$4} END {print s}' $1`;

A2=$[`wc -l < $2`-1];
B2=`awk '{s+=$4} END {print s}' $2`;

echo "Для файла $1:";
echo "$[(A1/B1)*100]%"
echo "Для файла $2:";
echo "$[(A2/B2)*100]%"
