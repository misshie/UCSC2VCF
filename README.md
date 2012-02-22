# UCSC2VCF
## WARNING!
The author knows that if a broken VCF file is used for an important project to filter-out potential variants, it results a CATASTROPHE. This is a MIT-licensed software. Use at your own risk ;-)

## Motivation
UCSC curates NCBI's dbSNP data before release at the UCSC Genome Database. However, only NCBI releases the dbSNP information in the VCF format. If you need UCSC-curated dbSNP information (dbSNP13x, dbSNP13xCommon, etc.) in the VCF format, you have to do that by yourself.

## Dependency
 gem install bio
 gem install bio-ucsc-api
 gem install striuct

## Prepare file
1. download reference sequence in the 2bit file format. For example, the human hg19 reference sequence is available at http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
2. download UCSC's dbSNP table using "Table Browser". Select "Group: Variation and Repeats" and "All SNPs (135)" or your favorite database. Then select "genome", select "all fields from selected table", and download (using gzip compress is better).
3. The VCF header is hard-coded in the script. Please change it if you do not use dbSNP135. 

## How to use
    ucsc2vcf -r reference.2bit --snv --indel > results.vcf

* `--snv` outputs "single nucleotide -> single nucleotide substitutions". Can be combined with `--indel`.
* `--indel` outputs " '-' to nucleotide(s) or nucleotide(s) to '-' substitutions". Can be combined with `--snv`.
* Other variations such as large indels (allels are not [-ATGC]) are all ignored.

## Copyright
 MISHIMA, Hiroyuki (hmishima at nagasaki-u.ac.jp, @mishimahryk at Twitter), 2012
