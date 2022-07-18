# MethylationCowparison

This script has been developed to compare the methylation patterns detected using Illumina's bisulfite sequencing and Oxford Nanopore Technologies methylation calling.
Input files need to be obtained using Bismarks for the BS and Nanopolish for ONT sequencing. 

Software is written in fortran90, and because of that the files need to be formatted as follows:

## Bisulfite file
Input file needs to be called Bisulfite.txt. 
File must not have a header.
5 columns with chromosome (integer), start site (integer), called sites (integer), methylated sites (integer), frequency (float up to 100).

Sexual chromosome need to be renamed as 30 and mitocondrial chromosome as 31. This can be done like:
> sed -r 's/X/30/g' InputFile.txt | sed -r 's/MT/31/g' | sort -n -k 1 > Bisulfite.txt
 
  
## Nanopore file
Input file needs to be called methylation_renamed.tsv
File must not have a header.
8 columns in the following order: chromosome (integer), start site (integer), end site (integer), CpGs in region (integer), called sites (integer), methylated sites (integer), methylation frequency (float up to 1), sequence (character).

Again, sexual chromosome need to be renamed as 30 and mitocondrial chromosome as 31. This can be done like:
> sed -r 's/X/30/g' InputFile.txt | sed -r 's/MT/31/g' | sort -n -k 1 > methylation_renamed.tsv


## Running
To compile this program:
> gfortran MethylationCowparison.f90 -o MethylationCowparison.x

To run it:
> MethylationCowparison.x

Methylation input files need to be in the same folder as the executable.


## Output file
Outputfile will be called ComparedSites.tsv and consists of 6 columns: 

- Chromosome in the Bisulfite file.
- Chromosome in the ONT file.
- Methylation frequency up to 1 in the BS file.
- Methylation frequency in the ONT file.
- Coverage of that cytosine in the BS file.
- Coverage of that cytosine in the ONT file. 

## Heatmap
(Adapted from Nanopolish)

> library(ggplot2)
> library(RColorBrewer)
> MethylationCowparison <- read.table("ComparedSites.tsv", header=T)
> colnames(MethylationCowparison) <- c("chrBS","chrONT","freqBS","freqONT","covBS","covONT")
> colours <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
> coloursr <- colours(32)
> 
> correlationSites <- cor(MethylationCowparison$freqBS, MethylationCowparison$freqONT)
> 
> title <- sprintf("Number of sites = %d r = %.3f", nrow(MethylationCowparison), correlationSites)
> ggplot(MethylationCowparison, aes(freqBS, freqONT)) +
>     geom_bin2d(bins=25) + scale_fill_gradientn(colors=coloursr, trans="log10") +
>     xlab("Bisulfite Methylation Frequency") +
>     ylab("Nanopolish Methylation Frequency") +
>     theme_bw(base_size=20) +
>     ggtitle(title)
