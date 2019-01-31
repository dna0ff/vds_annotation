#!/bin/bash

docker run -t -i \
  -v /data/d/vep_data:/opt/vep/.vep \
  -v /data/d/vep_data:/home/vep/.vep \
  -v /data/d/GRCh37:/GRCh37 \
  -v /data/d/pipe/_/tmp:/data \
  ensemblorg/ensembl-vep:release_90.10 \
  ./vep \
  --format vcf -i /data/05_mgrb_split_variants.vcf.bgz \
  --tab -o /data/05_mgrb_split_variants.vep.tab \
  --fasta /GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
  --no_stats --force_overwrite --assembly GRCh37 --minimal --everything --pick --gencode_basic --cache --offline --fork 56

