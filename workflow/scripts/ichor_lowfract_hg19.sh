ichor_lowfract() {
# Runs ichorCNA to generate tumor fraction
#  See https://doi.org/10.1038/s41467-017-00965-y
#
# Input parameters
#  $1 = input wig
#  $2 = output directory
#
# Steps
##
## Setup in-function parameters    
base=$(basename -s .wig $1)
##
## Check for inputs and outputs
if [ ! -f $1 ]; then
   echo "No input wig found"
elif [ $2/${base}.RData -nt $1 ]; then
   echo "wig for ${base} already processed in ichor"
else
   Rscript /opt/ichorCNA/scripts/runIchorCNA.R \
           --id $base \
           --WIG $1 \
           --gcWig /opt/ichorCNA/inst/extdata/gc_hg19_1000kb.wig \
           --normal "c(0.95, 0.99, 0.995, 0.999)" \
           --ploidy "c(2)" \
           --maxCN 3 \
           --estimateScPrevalence FALSE \
           --scStates "c()" \
           --outDir $2
fi
}
#
