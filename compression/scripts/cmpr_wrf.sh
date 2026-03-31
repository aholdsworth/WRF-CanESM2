#!/bin/bash

# helper script for compressing netcdf files

infile=$1
OUTPUT_DIR=$2

#OUTPUT_DIR="/gpfs/fs7/dfo/hpcmc/pfm/spfm000/DATA/WRF-CanESM2-Cmpr/"
# Check if files exist (in case glob doesn't match)
[ -e "$infile" ] || continue

# Extract the base filename
basename=$(basename "$infile")

# Create output filename (replace wrfout with CanESM2-WRF)
outfile="$OUTPUT_DIR/${basename/wrfout/CanESM2-WRF}.nc"


# Skip if output already exists
if [ -f "$outfile" ]; then
echo "[$current/$total] SKIP: $(basename "$outfile") already exists"
exit
fi

echo "[$current/$total] Processing: $basename"
echo "           Output: $(basename "$outfile")"

ncks -7 -L 4 \
  --ppc default=4 \
  --ppc Q2=6 --ppc QVAPOR=6 \
  --ppc PSFC=6 --ppc T2=5 --ppc T=5 --ppc THM=5 \
  --ppc U10=6 --ppc V10=6 --ppc U=6 --ppc V=6 \
  --ppc P=6 --ppc P_HYD=6 \
  --ppc PH=6 --ppc W=6 \
  --ppc HFX=5 --ppc LH=5 \
  --ppc RAINNC=5 --ppc SNOWNC=5 --ppc SNOW=5 \
  -x -v DX2D,CF1,CF2,CF3,RAINSH,HAILNC,ISEEDARRAY_SPP_CONV,ISEEDARRAY_SPP_LSM,ISEEDARRAY_SPP_PBL,ISEEDARR_RAND_PERTURB,ISEEDARR_SKEBS,ISEEDARR_SPPT,I_ACLWDNT,I_ACLWDNTC,RESM,ZETATOP,SSTSK,SWNORM,THIS_IS_AN_IDEAL_RUN,SAVE_TOPO_FROM_REAL \
  "$infile" "$outfile"

# Check if compression was successful
if [ $? -eq 0 ] && [ -f "$outfile" ]; then
# Print file sizes for comparison
insize=$(du -h "$infile" | cut -f1)
outsize=$(du -h "$outfile" | cut -f1)
echo "           SUCCESS: $insize -> $outsize"
else
echo "           FAILED"
# Remove failed output file if it exists
[ -f "$outfile" ] && rm "$outfile"
fi

echo ""
