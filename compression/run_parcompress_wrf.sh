datestring=`date +%Y%m%d-%H%M%S`
WRK_DIR="/gpfs/fs7/dfo/hpcmc/pfm/spfm000/amh001/WRF-Tests/"
INPUT_DIR="/home/spfm000/evg000/CanESM2_WRF_runs/historical_r1i1p1_1986/WRF/"
OUTPUT_DIR="/gpfs/fs7/dfo/hpcmc/pfm/spfm000/DATA/WRF-CanESM2-Cmpr/1986/"
FILELIST=${WRK_DIR}/${datestring}
(cd ${INPUT_DIR} && ls -1S wrfout_d03* > $FILELIST)
bash ${WRK_DIR}/parcompress_wrf.sh ${INPUT_DIR} ${OUTPUT_DIR} ${FILELIST}
