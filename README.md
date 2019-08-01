# genomon2vcf
Sample file is RK113_C01.genomon_mutation.result.filt.txt.

```
usage: genomon2vcf.py [-h] [--out VCF] [--ref FAS] [--sample ID] genomon_file

Convert Genomon2 SNV/INDEL calls to VCF.

positional arguments:
  genomon_file        Genomon2 SNV/INDEL file for a single sample

optional arguments:
  -h, --help          show this help message and exit
  --out VCF, -o VCF   Output VCF file [default: stdout]
  --ref FAS, -r FAS   FASTA file of reference human genome
  --sample ID, -s ID  Sample ID
```
