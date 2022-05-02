# MEGAdt
Design tools for Mutagenesis by Template-guided Amplicon Assembly (MEGA)

We have tested these scripts on Linux and MacOS

## Dependencies

* Python 3.9.7
	- os
	- sys
	- argparse
	- [Biopython](https://biopython.org/) - for melting temperature calculation

### Description
```
usage: MEGAdt.py [-f FASTA] [-i VARIANT] [-o OUTPUT] [Optional arguments]

Design tools for Mutagenesis by Template-guided Amplicon Assembly (MEGA). See details in https://github.com/hym0405/MEGAdt

Input and output [Required]:
  -f FASTA, --fasta FASTA
                        Input template sequences in FASTA format (See example: ./test_input/input_seq.fasta)
  -i VARIANT, --variant VARIANT
                        Input variants information in tab-separated table (See example: ./test_input/variant_info.tsv)
  -o OUTPUT, --output OUTPUT
                        Prefix of output files; [output].variant.fasta and [output].oligo.tsv will be generated

Oligo design [Optional]:
  -el END_LENGTH, --end_length END_LENGTH
                        Minimum number of perfectly matched bases (nt) at the end of mutagenesis oligo (both 3' and 5')' [default: 10]
  -gc3 {0,1}, --gc_clamp3 {0,1}
                        # of GC bases (nt) required at 3' end of mutagenesis oligo. The value should be either 1 or 0 [default: 1]
  -Tm5 MELTING5, --melting5 MELTING5
                        Minimum melting temperature (°C) at 5' end of mutagenesis oligos [default: 35]
  -Tm3h MELTING3_HIGH, --melting3_high MELTING3_HIGH
                        Upper limit of gradient melting temperature (°C) at 3' end of mutagenesis oligos [default: 65]
  -Tm3l MELTING3_LOW, --melting3_low MELTING3_LOW
                        Lower limit of gradient melting temperature (°C) at 3' end of mutagenesis oligos [default: 50]

Warning information [Optional]:
  -wg WARNING_GAP, --warning_gap WARNING_GAP
                        Minimum number of gaps (nt) between oligo to print warnings [default: 5]
  -wl WARNING_LENGTH, --warning_length WARNING_LENGTH
                        Minimum length of oligos to print warnings [default: 60]
  -h, --help            show this help message and exit
```

### Input format

**Input template sequences:** Input template sequences in FASTA format.

****[example: ./test_input/input_seq.fasta]****

```
>rsgA
CCATTGTTTTGTCGTTCCTGAT...
...
>pheS
ATGTCACATCTCGCAGAACTGG...
...
>AAVcap
ATGGCTGCCGATGGTTATCTTC...
...
```

**Input variants information:** Information of desired mutations of each variant. Variants can be SNPs, multiple-SNPs, insertions or deletions but at least one base is required in Reference_base (column-4) and Alternative_base (column-5).

****[Important] Positions of mutations (column-3) are 1-based****

****[Important] Template of variants (column-2) should match with input template sequences****

****[example: ./test_input/variant_info.tsv]****

```
Variant		Template	Position	Reference_base	Alternative_base
## single-site mutations
pheS_variant	pheS		397		A		G
pheS_variant	pheS		559		A		G
...
## multiple-sites mutations
rsgA_variant	rsgA		255		TAA		GAG
rsgA_variant	rsgA		408		CGA		GCT
...
## insertion and deletion
AAVcap_variant	AAVcap		416		C		CTGA
AAVcap_variant	AAVcap		570		GCCA		G
...
```

### Output format

****Results of mutagensis oligo designs ([output_prefix].oligo.tsv) and final sequences of desired variants ([output_prefix].variant.fasta)****

****[example of mutagensis oligo designs: ./test_output/output.oligo.tsv]****
```
Variant		Template	Primer			Start_site	Sequence				Length	Tm_5p	Tm_3p	Mutations	Notes				Warnings
rsgA_variant	rsgA		rsgA_variant-P0		1		CCATTGTTTTGTCGTTC			17	0	51.5			Extension F-primer	
rsgA_variant	rsgA		rsgA_variant-P1		244		TTACACAGACCgagATAGTCATGGAATTCGAC	32	35.4	55.3	255_TAA_GAG		
rsgA_variant	rsgA		rsgA_variant-P2		398		GACCGAGCCCgctGTTGTCAGAGATATCGTTG	32	43.8	57.3	408_CGA_GCT			
...
pheS_variant	pheS		pheS_variant-P5		326		TTACCCGTACCgTCGACCGTgTCGAAAGTTTCTTCGGTG	39	38.9	59.2	337_A_G;346_A_G	Contains multiple mutations	
pheS_variant	pheS		pheS_variant-P6		387		CGGGCCGGAAgTCGAAGACGATTATCATAACTTC	34	47.6	61	397_A_G		
...
AAVcap_variant	AAVcap		AAVcap_variant-P2	405		TGTTAAGACGGctgaTCCGGGAAAAAAGAGG		31	36.5	56.5	416_C_CTGA
AAVcap_variant	AAVcap		AAVcap_variant-P3	560		CTCTCGGACAgCCAGCAGCCCCCTC		25	36	61.4	570_GCCA_G
...
```
