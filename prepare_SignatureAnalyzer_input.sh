#!/bin/bash

bcftools=/usr/local/package/samtools/1.9/bin/bcftools
echo -e "Tumor_Sample_Barcode\tVariant_Type\tCHROM\tPOS\tReference_Allele\tTumor_Seq_Allele2\tVAF\tref_context"

for vcf in run-genomon2vcf.out/*.somatic.vcf.gz; do
    sample=`basename $vcf .somatic.vcf.gz`
    $bcftools query -f "${sample}\tSNP\t%CHROM\t%POS\t%REF\t%ALT\t%VAF\t%Context\n" $vcf | awk '$8 != "."'
done
