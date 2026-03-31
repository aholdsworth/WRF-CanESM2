#!/bin/bash

# Compress lots of files in parallel; requires gnu parallel and cmpr2.sh
PARALLEL=/home/stod000/bin/gnuparallel
COMPRESS=/home/spfm000/bin/cmpr_wrf.sh
set -e

# Paths must be full paths or this will probably not work
if [ "$#" -lt 3 ]; then
  echo "usage: $0 <DIR_TMP_CDF> <DIR_ARCHI_CDF> <FILELIST>"
  exit
fi
DIR_TMP_CDF=$1
DIR_ARCHI_CDF=$2
FILELIST=$3
echo $FILELIST
if [ ! -s $FILELIST ]; then
  echo "File list missing or empty"
  exit
fi

SCRIPT=${FILELIST}.sh
LOGFILE=${FILELIST}.log
TMPDIR=${DIR_TMP_CDF}/TMP

cd $DIR_TMP_CDF

sed \
-e s#DIR_TMP_CDF#$DIR_TMP_CDF#g \
-e s#DIR_ARCHI_CDF#$DIR_ARCHI_CDF#g \
-e s#FILELIST#$FILELIST#g \
-e s#TMPDIR#$TMPDIR#g \
-e s#LOGFILE#$LOGFILE#g \
<<EOF > $SCRIPT
#!/bin/bash
#SBATCH --export=USER,LOGNAME,HOME,MAIL,PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
#SBATCH --job-name=netcdf_parcompress
#SBATCH --output=LOGFILE
#SBATCH --qos=low
#SBATCH --account=dfo_pfm
#SBATCH --partition=standard
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=32000M
#SBATCH --comment="image=registry.maze.science.gc.ca/ssc-hpcs/generic-job:ubuntu22.04"

 
cd DIR_TMP_CDF
mkdir -p TMPDIR
$PARALLEL --tmpdir TMPDIR --delay 0.01 -j 24 -a FILELIST $COMPRESS {}
EOF

jobsub -c gpsc7 $SCRIPT

