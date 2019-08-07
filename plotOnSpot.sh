#!/bin/bash

## Script to plot the imaginary-time Green's functions. Terminal variables are: 1) filename and 2) step number for plotting (step along beta).

filename="$1"
steps=$2

if [ -z $filename ]; then
        echo "no file input!"
        exit 1
fi

#arraynums=$(echo $filename | sed -e 's/[^0-9]*//g')
#arraynums=$(echo $filename | sed -e 's/\(.*[a-zA-Z]_\)\([0-9]*\)/\2/g')

# The format of the string has to be: .*beta_betainit_betamax_betainc.dat
betainit=$(echo $filename | sed -e 's/\(.*[beta]_\)\([0-9]*\)/\2/g' | sed -e 's/\..*//g' | awk -F "_" '{print $1}')
betamax=$(echo $filename | sed -e 's/\(.*[beta]_\)\([0-9]*\)/\2/g' | sed -e 's/\..*//g' | awk -F "_" '{print $2}')
betainc=$(echo $filename | sed -e 's/\(.*[beta]_\)\([0-9]*\)/\2/g' | sed -e 's/\..*//g' | awk -F "_" '{print $3}')
echo -e "$betainit\n$betamax\n filename: $filename"

numoflines=$( cat $filename | wc -l )

## Function to skip some lines in plotting.
skipfile(){
        stringtmp="\"< awk '(NR==$1){print;}' $2\""
        echo $stringtmp
}

echo -e "Number of lines $numoflines\n"

for l in $(seq 1 $steps $numoflines)
do
       jump=$( skipfile $l $filename )
       echo $jump
       if [ $l -eq 1 ]; then
        	gnuplot -e "set macros; set title \"G(tau)\"; plot $jump matrix w lp; pause mouse key"
       else
	gnuplot -e "set macros; set title \"G(tau)\"; plot $jump matrix w lp; pause mouse key"
	#continue
       fi
done
