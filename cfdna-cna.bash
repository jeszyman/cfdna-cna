cd ~/repos/cfdna-cna
ls
mkdir -p tests/full-project/analysis/cna-wig

cd ~/repos/cfdna-cna
git clone git@github.com:broadinstitute/ichorCNA ichorCNA.v0.2.0
cd ichorCNA.v0.2.0
git fetch --depth 1 origin tag v0.2.0
git checkout tags/v0.2.0

cd ~/repos/cfdna-cna
conda activate base
conda create --name ext-cfdna-cna -y

curl -o /tmp/basecamp_env.yaml https://raw.githubusercontent.com/jeszyman/basecamp/v1.0.0/basecamp_env.yaml && mamba env update --name ext-cfdna-cna --file /tmp/basecamp_env.yaml

curl -o /tmp/biotools_env.yaml https://raw.githubusercontent.com/jeszyman/biotools/v1.0.0/biotools_env.yaml && mamba env update --name ext-cfdna-cna --file /tmp/biotools_env.yaml

mamba env update --name ext-cfdna-cna --file config/ext-cfdna-cna.yaml

singularity shell --bind /mnt ~/sing_containers/cfdna_wgs.1.0.0.sif

# Clear bam directory if present
if [ -r test/bam ]; then \rm -rf test/bam; fi
mkdir -p test/bam

# Create small bam files to store in repo. Subsample real bams to ~100 Mb.
sambamba view -s .005 -f bam -t 36 /mnt/ris/aadel/mpnst/bam/lib070_dedup_sorted.bam > test/bam/lib001_hg19.bam
sambamba view -s .005 -f bam -t 36 /mnt/ris/aadel/mpnst/bam/lib054_dedup_sorted.bam > test/bam/lib002_hg19.bam
sambamba view -s .005 -f bam -t 36 /mnt/ris/aadel/mpnst/test/bam/new_HiSeq15_L002001_ACAC_extract_ds20.bam > test/bam/lib003_hg38.bam
sambamba view -s .005 -f bam -t 36 /mnt/ris/aadel/mpnst/test/bam/new_HiSeq15_L002001_ATCG_extract_ds20.bam > test/bam/lib004_hg38.bam

sambamba view -s 0.01 -f bam -t 4 /mnt/ris/aadel/mpnst/bam/cfdna_wgs/ds/lib105_ds10.bam > test/bam/lib005.bam
sambamba view -s 0.01 -f bam -t 4 /mnt/ris/aadel/mpnst/bam/cfdna_wgs/ds/lib205_ds10.bam > test/bam/lib006.bam



for file in test/bam/*.bam; do samtools index $file; done

cd ~/repos/cfdna-cna
conda activate ext-cfdna-cna

Rscript ./ichorCNA.v0.2.0/scripts/runIchorCNA.R \
         --id TEST \
         --WIG ./tests/wigs/test.wig \
         --gcWig ./ichorDNA.v0.2.0/inst/extdata/gc_hg38_1000kb.wig \
         --mapWig ./ichorDNA.v0.2.0/inst/extdata/map_hg38_1000kb.wig \
         --centromere ./ichorCNA.v0.2.0/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
         --normal "c(0.95, 0.99, 0.995, 0.999)" \
         --normalPanel ./ichorCNA.v0.2.0/inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds \
         --ploidy "c(2)" \
         --maxCN 3 \
         --estimateScPrevalence FALSE \
         --scStates "c()" \
         --outDir ./test/ichor

  # mkdir /tmp/ichor_out
  # singularity shell ~/sing_containers/mpnst.sif

  Rscript /opt/ichorCNA/scripts/runIchorCNA.R --id tumor_sample \
  --WIG /tmp/test.wig --ploidy "c(2,3)" --normal "c(0.5,0.6,0.7,0.8,0.9)" --maxCN 5 \
  --gcWig /opt/ichorCNA/inst/extdata/gc_hg38_1000kb.wig \



  --includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" \
  --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
  --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --outDir /tmp/ichor_out

  #mkdir -p /tmp/ichor_out
  #singularity shell ~/sing_containers/mpnst.sif

  # Notes
  ##
  ## Will overwrite target files with a warning
  ##
  ##


  Rscript /opt/ichorCNA/scripts/runIchorCNA.R --id tumor_sample \
    --WIG ~/repos/cfdna-cna/test/wig/lib002_hg19_frag90_150.wig --ploidy "c(2,3)" --normal "c(0.5,0.6,0.7,0.8,0.9)" --maxCN 5 \
    --gcWig /opt/ichorCNA/inst/extdata/gc_hg19_1000kb.wig \
    --mapWig /opt/ichorCNA/inst/extdata/map_hg19_1000kb.wig \
    --centromere /opt/ichorCNA/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt \
    --normalPanel /opt/ichorCNA/inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds \
    --includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" \
    --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
    --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --outDir /tmp/ichor_out

  # mkdir /tmp/ichor_out
  # singularity shell ~/sing_containers/mpnst.sif

  Rscript ./workflow/scripts/MOD_runIchorCNA.R --id tumor_sample \
    --WIG ~/repos/cfdna-cna/test/wig/lib002_frag90_150.wig --ploidy "c(2,3)" --normal "c(0.5,0.6,0.7,0.8,0.9)" --maxCN 5 \
    --gcWig /opt/ichorCNA/inst/extdata/gc_hg19_1000kb.wig \
    --mapWig /opt/ichorCNA/inst/extdata/map_hg19_1000kb.wig \
    --centromere /opt/ichorCNA/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt \
    --normalPanel /opt/ichorCNA/inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds \
    --includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" \
    --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
    --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --outDir /tmp/ichor_out

      Rscript /opt/ichorCNA/scripts/runIchorCNA.R --id tumor_sample \
        --WIG /tmp/test_hg19.wig --ploidy "c(2,3)" --normal "c(0.5,0.6,0.7,0.8,0.9)" --maxCN 5 \
        --gcWig /opt/ichorCNA/inst/extdata/gc_hg19_1000kb.wig \
        --mapWig /opt/ichorCNA/inst/extdata/map_hg19_1000kb.wig \
        --centromere /opt/ichorCNA/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt \
        --normalPanel /opt/ichorCNA/inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds \
        --includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" \
        --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
        --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --outDir /tmp/ichor_out_test

  Rscript workflow/scripts/MOD_runIchorCNA.R --id lib003_hg38_frag90_150 \
          --WIG test/wig/lib003_hg38_frag90_150.wig \
          --gcWig /opt/ichorCNA/inst/extdata/gc_hg38_1000kb.wig \
          --mapWig /opt/ichorCNA/inst/extdata/map_hg38_1000kb.wig \
          --centromere /opt/ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt \
          --normal "c(0.95, 0.99, 0.995, 0.999)" \
          --ploidy "c(2)" \
          --maxCN 3 \
          --estimateScPrevalence FALSE \
          --scStates "c()" \
          --outDir test/ichor_hg38



          # --normalPanel /opt/ichorCNA/inst/extdata/HD_ULP_PoN_hg38_1Mb_median_normAutosome_median.rds \
  #                                                   HD_ULP_PoN_hg38_1Mb_median_normAutosome_median.rds

#singularity shell --bind /mnt:/mnt ~/sing_containers/mpnst.sif


/opt/hmmcopy_utils/bin/readCounter \
 --window 1000000 \
 --quality 20 \
 --chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX" \
 ~/repos/cfdna-cna/test/bam/lib004_hg38.bam > /tmp/test.wig
