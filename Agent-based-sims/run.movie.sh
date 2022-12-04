#!/bin/bash

# geometric parameters
alpha0=1.0
alphamax=2.0
Lx=32.0
att=0.0
s=0.0

# rates
rateWT=1d0
b=4e3

# steps
steps_burn=400000
steps_decorr=40000
trials=10

# skip steps
layerskip=25
dataskip=10000
prodskip=2000
restskip=10000
dt=2e-5

# growth layer widths
layerwidth=1.0
layerdepth=13.0
propdepth=5.0
bounddepth=3.0

# run parameters
desync=0.4
seed0=-9

# output files
file1=outs/prod.movie.dat
file2=restart_WT.movie.dat
file3=restart_MUT.movie.dat
file4=mut_sizes.movie.dat
file5=out_config.movie.dat
file7=out.movie.dat
file8=width.movie.dat
file9=front_line.movie.dat

# logical variables
movie=.true.
restart=.true.
bottom=.false.

# output directory
outdir=outs

./budding.damped.linear_front.serial_mut.sh \
        -a $alpha0 \
	-b $alphamax \
        -c $Lx \
        -d $att \
        -e $s \
        -f $rateWT \
        -g $b \
        -h $steps_burn \
        -i $steps_decorr \
        -j $trials \
        -k $layerskip \
        -l $dataskip \
        -m $prodskip \
        -n $restskip \
        -o $dt \
        -p $layerwidth \
        -q $layerdepth \
        -r $propdepth \
        -s $bounddepth \
        -t $desync \
        -u $seed0 \
        -v $file1 \
        -w $file2 \
        -x $file3 \
        -y $file4 \
        -z $file5 \
        -A $file7 \
        -B $file8 \
        -C $file9 \
        -D $movie \
        -E $restart \
        -F $bottom \
        -G $outdir

