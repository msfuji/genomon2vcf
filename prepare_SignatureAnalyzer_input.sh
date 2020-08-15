#!/bin/bash

SAMTOOLS_DIR=/usr/local/package/samtools/1.9/bin
echo -e "Tumor_Sample_Barcode\tVariant_Type\tCHROM\tPOS\tReference_Allele\tTumor_Seq_Allele2\tVAF\tref_context"

for vcf in *.vcf; do
    sample=`basename $vcf .vcf`
    $SAMTOOLS_DIR/bgzip $vcf
    $SAMTOOLS_DIR/tabix -p vcf ${vcf}.gz
    $SAMTOOLS_DIR/bcftools query -f "${sample}\tSNP\t%CHROM\t%POS\t%REF\t%ALT\t%VAF\t%Context\n" ${vcf}.gz | awk '$8 != "."'
done
