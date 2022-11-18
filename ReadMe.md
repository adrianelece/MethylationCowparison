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
 
Remove scaffolds: 
> sed -i '/^NW/d' Bisulfite.txt

## Nanopore file
Input file needs to be called methylation_renamed.tsv
File must not have a header.
8 columns in the following order: chromosome (integer), start site (integer), end site (integer), CpGs in region (integer), called sites (integer), methylated sites (integer), methylation frequency (float up to 1), sequence (character).

Remove scaffolds: 
> sed -i '/^NW/d' Bisulfite.txt

Again, sexual chromosome need to be renamed as 30 and mitocondrial chromosome as 31. This can be done like:
> sed -r 's/X/30/g' InputFile.txt | sed -r 's/MT/31/g' | sort -n -k 1 > methylation_renamed.tsv


## Running
To compile this program:
> gfortran MethylationCowparison.f90 -o MethylationCowparison.x

To run it:
> ./MethylationCowparison.x

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
```sh
library(ggplot2)
library(RColorBrewer)
setwd("dir/")
#library(scales)
#library(peRReo)
#paleta<-latin_palette("rosalia",n=100)

MethylationCowparison <- read.table("ComparedSites.tsv", header=F)
colnames(MethylationCowparison) <- c("chrBS","chrONT","freqBS","freqONT","covBS","covONT")
colours <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
coloursr <- colours(32)

samplen <- "X"
mincovONT <- 25
mincovBS <- 25


MethylationCowparisonfiltered <- MethylationCowparison[which(MethylationCowparison$covONT>=mincovONT&MethylationCowparison$covBS>=mincovBS),]
filtered <- which(MethylationCowparisonfiltered$freqONT==0|MethylationCowparisonfiltered$freqONT==1)

MethylationCowparisonfiltered2 <- MethylationCowparisonfiltered[-filtered,]
correlationSites <- cor(MethylationCowparisonfiltered2$freqBS, MethylationCowparisonfiltered2$freqONT)

title <- sprintf("Number of sites = %d r = %.3f sample = %s", nrow(MethylationCowparisonfiltered2), correlationSites,samplen)
ggplot(MethylationCowparisonfiltered2, aes(x=freqBS, y=freqONT)) +
  geom_bin2d(bins=25) + scale_fill_gradientn(colors=coloursr, trans="log10") +
# geom_bin2d(bins=15) + scale_fill_gradientn(colors=paleta, trans="log10",limits=c(50,15000),oob=censor,na.value = "white") +
  xlab("Bisulfite Methylation Frequency") +
  ylab("Nanopolish Methylation Frequency") +
  theme_bw(base_size=20) +
  ggtitle(title)

```



# Weighted genome wide methylation

## To transform the files from modbam2bed into the necessary input format:

```
awk '($11!="nan"&&$11!=0){print $1"\t"$2"\t"$6"\t""CG""\t"$11*$10/100"\t"$10}' X.freq.txt > 5014_1x_physalia.tsv
```


```
#!/usr/bin/env Rscript

# load packages
library(data.table)
library(ggplot2)

# set working directory
DirectoryPath <- "C:/Users/lop025/OneDrive - CSIRO/Documents/INIA/mod2bam/"
setwd(DirectoryPath)

# chromosome	position	strand(+/-)	mc_class(CG/CHG/CHH)	methylated_bases(number_methylated_C_reads)	total_bases(total_number_reads_mapping)


###################################################################################################
#                   Plotting Gene methylation metaplot
###################################################################################################
# loading data
gff3 <- fread("genesBtau.gff3", h=T)
data1 <- fread("3130_1x_physalia.tsv", h=T)
data1$methylated_bases = round(data1$methylated_bases)
data5 = data1[which(data1$total_bases>5),]
data7 = data1[which(data1$total_bases>7),]
data10 = data1[which(data1$total_bases>10),]


data=data1
# loop on all genes
for (i in seq(1, nrow(gff3))) {
  # gene info
  CurrentGene <- gff3$GeneName[i]
  CurrentChr <- gff3$chr[i]
  Currentmin <- gff3$start[i] - 4001
  Currentmax <- gff3$stop[i] + 4001
  Currentorientation <- gff3$strand[i]
  
  # extract gene methylation and check there is something
  GeneMethylation <- data[(chrom==CurrentChr)&(pos>=Currentmin)&(pos<=Currentmax)]
  if (nrow(GeneMethylation)>0) {
    # if gene is on minus strand, rewrite positions in increasing order for practical calculations
    GeneMethylation[, Newpos := if(Currentorientation=='+'){pos}else{Currentmax-pos+Currentmin}, by=1:nrow(GeneMethylation)]
    
    # define window boundaries
    geneWindow <- floor((Currentmax - Currentmin - 8002 + 1) / 20)
    breaks <- c(seq(from=Currentmin, to=Currentmin+3800, by=200), seq(from=Currentmin+4001, to=Currentmin+4001+19*geneWindow, by=geneWindow), Currentmax-4001, seq(from=Currentmax-3800, to=Currentmax, by=200))
    breaks[41] <- Currentmax-4000
    
    # loop on each window
    for (w in seq(1, length(breaks)-1)) {
      # extract window methylation info
      WindowMethylation <- GeneMethylation[(Newpos>=breaks[w])&(Newpos<breaks[w+1])]
      
      if (nrow(WindowMethylation)>0) {
        summary <- WindowMethylation[, .(SumMethylatedReads = sum(methylated_bases), TotalReads = sum(total_bases)), by=mc_class]
        if (nrow(summary[mc_class == "CG"]) == 0) {
          summary <- as.data.table(rbind(summary, list("CG", 0, 0)))
        }
        summary[, window := w]
        summary[, GeneName := CurrentGene]
        
        # add gene info to output file
        fwrite(summary, file="Gene_windows1.txt", sep="\t", append=TRUE)
      }
    }
  }
}


# load windows data
Windows1 <- fread("Gene_windows1.txt", h=T)
Windows1
colnames(Windows1)= c("mc_class",  "SumMethylatedReads", "TotalReads", "window", "GeneName")




######

data=data5
# loop on all genes
for (i in seq(1, nrow(gff3))) {
  # gene info
  CurrentGene <- gff3$GeneName[i]
  CurrentChr <- gff3$chr[i]
  Currentmin <- gff3$start[i] - 4001
  Currentmax <- gff3$stop[i] + 4001
  Currentorientation <- gff3$strand[i]
  
  # extract gene methylation and check there is something
  GeneMethylation <- data[(chrom==CurrentChr)&(pos>=Currentmin)&(pos<=Currentmax)]
  if (nrow(GeneMethylation)>0) {
    # if gene is on minus strand, rewrite positions in increasing order for practical calculations
    GeneMethylation[, Newpos := if(Currentorientation=='+'){pos}else{Currentmax-pos+Currentmin}, by=1:nrow(GeneMethylation)]
    
    # define window boundaries
    geneWindow <- floor((Currentmax - Currentmin - 8002 + 1) / 20)
    breaks <- c(seq(from=Currentmin, to=Currentmin+3800, by=200), seq(from=Currentmin+4001, to=Currentmin+4001+19*geneWindow, by=geneWindow), Currentmax-4001, seq(from=Currentmax-3800, to=Currentmax, by=200))
    breaks[41] <- Currentmax-4000
    
    # loop on each window
    for (w in seq(1, length(breaks)-1)) {
      # extract window methylation info
      WindowMethylation <- GeneMethylation[(Newpos>=breaks[w])&(Newpos<breaks[w+1])]
      
      if (nrow(WindowMethylation)>0) {
        summary <- WindowMethylation[, .(SumMethylatedReads = sum(methylated_bases), TotalReads = sum(total_bases)), by=mc_class]
        if (nrow(summary[mc_class == "CG"]) == 0) {
          summary <- as.data.table(rbind(summary, list("CG", 0, 0)))
        }
        summary[, window := w]
        summary[, GeneName := CurrentGene]
        
        # add gene info to output file
        fwrite(summary, file="Gene_windows5.txt", sep="\t", append=TRUE)
      }
    }
  }
}


# load windows data
Windows5 <- fread("Gene_windows5.txt", h=T)
Windows5
colnames(Windows5)= c("mc_class",  "SumMethylatedReads", "TotalReads", "window", "GeneName")


#####

data=data7
# loop on all genes
for (i in seq(1, nrow(gff3))) {
  # gene info
  CurrentGene <- gff3$GeneName[i]
  CurrentChr <- gff3$chr[i]
  Currentmin <- gff3$start[i] - 4001
  Currentmax <- gff3$stop[i] + 4001
  Currentorientation <- gff3$strand[i]
  
  # extract gene methylation and check there is something
  GeneMethylation <- data[(chrom==CurrentChr)&(pos>=Currentmin)&(pos<=Currentmax)]
  if (nrow(GeneMethylation)>0) {
    # if gene is on minus strand, rewrite positions in increasing order for practical calculations
    GeneMethylation[, Newpos := if(Currentorientation=='+'){pos}else{Currentmax-pos+Currentmin}, by=1:nrow(GeneMethylation)]
    
    # define window boundaries
    geneWindow <- floor((Currentmax - Currentmin - 8002 + 1) / 20)
    breaks <- c(seq(from=Currentmin, to=Currentmin+3800, by=200), seq(from=Currentmin+4001, to=Currentmin+4001+19*geneWindow, by=geneWindow), Currentmax-4001, seq(from=Currentmax-3800, to=Currentmax, by=200))
    breaks[41] <- Currentmax-4000
    
    # loop on each window
    for (w in seq(1, length(breaks)-1)) {
      # extract window methylation info
      WindowMethylation <- GeneMethylation[(Newpos>=breaks[w])&(Newpos<breaks[w+1])]
      
      if (nrow(WindowMethylation)>0) {
        summary <- WindowMethylation[, .(SumMethylatedReads = sum(methylated_bases), TotalReads = sum(total_bases)), by=mc_class]
        if (nrow(summary[mc_class == "CG"]) == 0) {
          summary <- as.data.table(rbind(summary, list("CG", 0, 0)))
        }
        summary[, window := w]
        summary[, GeneName := CurrentGene]
        
        # add gene info to output file
        fwrite(summary, file="Gene_windows7.txt", sep="\t", append=TRUE)
      }
    }
  }
}


# load windows data
Windows7 <- fread("Gene_windows7.txt", h=T)
Windows7
colnames(Windows7)= c("mc_class",  "SumMethylatedReads", "TotalReads", "window", "GeneName")


##############


data=data10
# loop on all genes
for (i in seq(1, nrow(gff3))) {
  # gene info
  CurrentGene <- gff3$GeneName[i]
  CurrentChr <- gff3$chr[i]
  Currentmin <- gff3$start[i] - 4001
  Currentmax <- gff3$stop[i] + 4001
  Currentorientation <- gff3$strand[i]
  
  # extract gene methylation and check there is something
  GeneMethylation <- data[(chrom==CurrentChr)&(pos>=Currentmin)&(pos<=Currentmax)]
  if (nrow(GeneMethylation)>0) {
    # if gene is on minus strand, rewrite positions in increasing order for practical calculations
    GeneMethylation[, Newpos := if(Currentorientation=='+'){pos}else{Currentmax-pos+Currentmin}, by=1:nrow(GeneMethylation)]
    
    # define window boundaries
    geneWindow <- floor((Currentmax - Currentmin - 8002 + 1) / 20)
    breaks <- c(seq(from=Currentmin, to=Currentmin+3800, by=200), seq(from=Currentmin+4001, to=Currentmin+4001+19*geneWindow, by=geneWindow), Currentmax-4001, seq(from=Currentmax-3800, to=Currentmax, by=200))
    breaks[41] <- Currentmax-4000
    
    # loop on each window
    for (w in seq(1, length(breaks)-1)) {
      # extract window methylation info
      WindowMethylation <- GeneMethylation[(Newpos>=breaks[w])&(Newpos<breaks[w+1])]
      
      if (nrow(WindowMethylation)>0) {
        summary <- WindowMethylation[, .(SumMethylatedReads = sum(methylated_bases), TotalReads = sum(total_bases)), by=mc_class]
        if (nrow(summary[mc_class == "CG"]) == 0) {
          summary <- as.data.table(rbind(summary, list("CG", 0, 0)))
        }
        summary[, window := w]
        summary[, GeneName := CurrentGene]
        
        # add gene info to output file
        fwrite(summary, file="Gene_windows10.txt", sep="\t", append=TRUE)
      }
    }
  }
}


# load windows data
Windows10 <- fread("Gene_windows10.txt", h=T)
Windows10
colnames(Windows10)= c("mc_class",  "SumMethylatedReads", "TotalReads", "window", "GeneName")

Windows1$Depth = 1
Windows5$Depth = 5
Windows7$Depth = 7
Windows10$Depth = 10

allcov = rbind(Windows1,Windows5,Windows7,Windows10)


# Summary over all genes
allcov$TotalReads = round(allcov$TotalReads)
allcov$SumMethylatedReads = round(allcov$SumMethylatedReads)
SummaryWindows <- allcov[TotalReads > 0, .(sum(SumMethylatedReads)/sum(TotalReads), binom.test(sum(SumMethylatedReads), sum(TotalReads))[[4]][1], binom.test(sum(SumMethylatedReads), sum(TotalReads))[[4]][2]), by=c("window", "mc_class","Depth")]

names(SummaryWindows) <- c("window", "mc_class", "Depth", "Average_Methylation", "CI_low", "CI_high")
SummaryWindows$Depth = as.character(SummaryWindows$Depth)

# plot metaplot

ggplot(data=SummaryWindows, aes(x = window, y = Average_Methylation, col=Depth)) + 
  theme_bw(base_size=20) +
  geom_point(aes(x = window, y = Average_Methylation, col=Depth), size = 1) +  
  geom_line(aes(x = window, y = Average_Methylation, col=Depth), size=0.5) +
  geom_ribbon(aes(ymin=CI_low, ymax=CI_high, col=Depth, fill=Depth), linetype=1, size=0.2, alpha=0.1, show.legend = FALSE) +
  labs(x=paste("Position", sep=''), y="Weighted methylation") + 
  guides(col=guide_legend(title="Coverage threshold")) +
  scale_x_continuous(breaks=c(1, 21, 40, 60), labels=c('-4000', 'TSS', 'TTS', '+4000'))

```
