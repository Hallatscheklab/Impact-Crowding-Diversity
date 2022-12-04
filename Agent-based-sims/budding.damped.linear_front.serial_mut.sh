#!/bin/bash

set -o errexit


# usage statement
usage() {
    cat << EOF

usage: $0 parameters

parameters:
  -a      Ratio of cell length to width at birth (alpha0)
  -b      Ratio of cell length to width at division (alphamax)
  -c      Length of box (Lx)
  -d      Range of attractive force (att)
  -e      Growth rate difference of mutant to WT cells (s)
  -f      Growth rate of WT cells (rateWT)
  -g      Friction coefficient for overdamped dynamics (b)
  -h      # steps to initialize before starting mutations (steps_burn)
  -i      # steps in between mutations (steps_decorr)
  -j      # mutations to introduce serially (trials)
  -k      # steps in between calulating distance to front (layerskip)
  -l      # steps in between outputting data (dataskip)
  -m      # steps in between outputting production file (prodskip)
  -n      # steps in between outputting restart file (restskip)
  -o      Time step (dt)
  -p      Numerical parameter for finding distance to front (layerwidth)
  -q      Depth of region where cells actively grow (layerdepth)
  -r      Depth of region where cells are pushed by forces (propdepth)
  -s      Depth of region where cell positions are fixed (bounddepth)
  -t      Noise on growth rate to desynronize cell cycles (desync)
  -u      Random seed (seed0)
  -v      Movie file (file1)
  -w      Restart file for burned config (file2)
  -x      Restart file for mutants (file3)
  -y      Mutant statistics (file4)
  -z      Mutant config (file5)
  -A      # cells in various layers (file7)
  -B      Width of mutant clone in growth layer (file8)
  -C      Binned front-line of colony (file9)
  -D      Logical parameter to output movie  (movie)
  -E      Logical parameter to initialize with restart file (restart)
  -F      Logical parameter to include bottom of colony (bottom)
  -G      Directory to output to (not used in Fortran program)
EOF
}


# make sure we have correct number of options
if [ $# -ne 66 ]
then
  echo "Invalid number of parameters"
  usage
  exit 1
fi


# read options
while getopts a:b:c:d:e:f:g:h:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:A:B:C:D:E:F:G: option
do
  case $option in
    a)
      alpha0=$OPTARG
      ;;
    b)
      alphamax=$OPTARG
      ;;
    c)
      Lx=$OPTARG
      ;;
    d)
      att=$OPTARG
      ;;
    e)
      s=$OPTARG
      ;;
    f)
      rateWT=$OPTARG
      ;;
    g)
      b=$OPTARG
      ;;
    h)
      steps_burn=$OPTARG
      ;;
    i)
      steps_decorr=$OPTARG
      ;;
    j)
      trials=$OPTARG
      ;;
    k)
      layerskip=$OPTARG
      ;;
    l)
      dataskip=$OPTARG
      ;;
    m)
      prodskip=$OPTARG
      ;;
    n)
      restskip=$OPTARG
      ;;
    o)
      dt=$OPTARG
      ;;
    p)
      layerwidth=$OPTARG
      ;;
    q)
      layerdepth=$OPTARG
      ;;
    r)
      propdepth=$OPTARG
      ;;
    s)
      bounddepth=$OPTARG
      ;;
    t)
      desync=$OPTARG
      ;;
    u)
      seed0=$OPTARG
      ;;
    v)
      file1=$OPTARG
      ;;
    w)
      file2=$OPTARG
      ;;
    x)
      file3=$OPTARG
      ;;
    y)
      file4=$OPTARG
      ;;
    z)
      file5=$OPTARG
      ;;
    A)
      file7=$OPTARG
      ;;
    B)
      file8=$OPTARG
      ;;
    C)
      file9=$OPTARG
      ;;
    D)
      movie=$OPTARG
      ;;
    E)
      restart=$OPTARG
      ;;
    F)
      bottom=$OPTARG
      ;;
    G)
      outdir=$OPTARG
      ;;
    *) 
      usage
      exit 1
      ;;
  esac
done

# change to output directory
rundir=$(pwd)
cd $outdir

# run program
time $rundir/budding.damped.linear_front.serial_mut.o <<EOF
  $alpha0
  $alphamax
  $Lx
  $att
  $s
  $rateWT
  $b
  $steps_burn
  $steps_decorr
  $trials
  $layerskip
  $dataskip
  $prodskip
  $restskip
  $dt
  $layerwidth
  $layerdepth
  $propdepth
  $bounddepth
  $desync
  $seed0
  $file1
  $file2
  $file3
  $file4
  $file5
  $file7
  $file8
  $file9
  $movie
  $restart
  $bottom
EOF
