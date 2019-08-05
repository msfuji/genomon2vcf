# genomon2vcf
Requires Python3. Sample file is RK113_C01.genomon_mutation.result.filt.txt.

```
usage: genomon2vcf.py [-h] [--out VCF] [--sample ID] genomon_file ref_genome

Convert Genomon2 SNV/INDEL calls to VCF.

positional arguments:
  genomon_file        Genomon2 SNV/INDEL file for a single sample
  ref_genome          FASTA file of reference human genome

optional arguments:
  -h, --help          show this help message and exit
  --out VCF, -o VCF   Output VCF file [default: stdout]
  --sample ID, -s ID  Sample ID
```
