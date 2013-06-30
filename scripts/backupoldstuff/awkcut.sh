#!/bin/bash

if [ $# -eq 1 ]
then
	fq=$1
	echo "path to file one and [ENTER]"
	read inputone
	echo "path to file two and [ENTER]"
	read inputtwo
	echo "path to output and [ENTER]"
	read output
elif [ $# -eq 3 ]
then 
	fq=400000
	inputone=$1
	inputtwo=$2
	output=$3
elif [ $# -eq 4 ]
then
	fq=$1
	inputone=$2
	inputtwo=$3
	output=$4
else [ $# -ne 4 ]
	echo "Sampling freq?"
	read fq
	echo "path to file one and [ENTER]"
	read inputone
	echo "path to file two and [ENTER]"
	read inputtwo
	echo "path to output and [ENTER]"
	read output

fi


awk -v frq=$fq 'NR>7 {print (NR-7)/frq " " $3 " " $4 ;}' $inputone > tmp1
awk 'NR>7 {print $3 " " $4 ;}' $inputtwo > tmp2

paste tmp1 tmp2 > tmp3
awk 'gsub("\r","",$0);' tmp3 > $output
rm tmp1 tmp2 tmp3

