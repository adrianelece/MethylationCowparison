# MethylationCowparison

This script has been developed to compare the methylation patterns detected using Illumina's bisulfite sequencing and Oxford Nanopore Technologies methylation calling.
Input files need to be obtained using Bismarks for the BS and Nanopolish for ONT sequencing. 

Software is written in fortran90, and because of that the files need to be formatted as follows:

## Bisulfite file
Input file needs to be called Bisulfite.txt
5 columns with chromosome (integer), start site (integer), called sites (integer), methylated sites (integer), frequency (float up to 100).

Sexual chromosome need to be renamed as 30 and mitocondrial chromosome as 31. This can be done like:
> sed -r 's/X/30/g' InputFile.txt | sed -r 's/MT/31/g' > Bisulfite.txt
 
  
## Nanopore file
Input file needs to be called methylation_renamed.tsv
8 columns in the following order: chromosome (integer), start site (integer), end site (integer), CpGs in region (integer), called sites (integer), methylated sites (integer), methylation frequency (float up to 1), sequence (character).

Again, sexual chromosome need to be renamed as 30 and mitocondrial chromosome as 31. This can be done like:
> sed -r 's/X/30/g' InputFile.txt | sed -r 's/MT/31/g' | sort -n -k 1 > methylation_renamed.tsv
