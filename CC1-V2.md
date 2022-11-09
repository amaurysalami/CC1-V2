TUTO
================

``` r
library(dada2)
```

    ## Loading required package: Rcpp

``` r
library(phyloseq)
library(DECIPHER)
```

    ## Loading required package: Biostrings

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     distance

    ## Loading required package: XVector

    ## Loading required package: GenomeInfoDb

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    ## Loading required package: RSQLite

    ## Loading required package: parallel

``` r
library(phangorn)
```

    ## Loading required package: ape

    ## 
    ## Attaching package: 'ape'

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     complement

``` r
library(ggplot2)
library(gridExtra)
```

    ## 
    ## Attaching package: 'gridExtra'

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     combine

``` r
miseq_path <- "/home/rstudio/MiSeq_SOP"
list.files(miseq_path)
```

    ##  [1] "F3D0_S188_L001_R1_001.fastq"   "F3D0_S188_L001_R2_001.fastq"  
    ##  [3] "F3D1_S189_L001_R1_001.fastq"   "F3D1_S189_L001_R2_001.fastq"  
    ##  [5] "F3D141_S207_L001_R1_001.fastq" "F3D141_S207_L001_R2_001.fastq"
    ##  [7] "F3D142_S208_L001_R1_001.fastq" "F3D142_S208_L001_R2_001.fastq"
    ##  [9] "F3D143_S209_L001_R1_001.fastq" "F3D143_S209_L001_R2_001.fastq"
    ## [11] "F3D144_S210_L001_R1_001.fastq" "F3D144_S210_L001_R2_001.fastq"
    ## [13] "F3D145_S211_L001_R1_001.fastq" "F3D145_S211_L001_R2_001.fastq"
    ## [15] "F3D146_S212_L001_R1_001.fastq" "F3D146_S212_L001_R2_001.fastq"
    ## [17] "F3D147_S213_L001_R1_001.fastq" "F3D147_S213_L001_R2_001.fastq"
    ## [19] "F3D148_S214_L001_R1_001.fastq" "F3D148_S214_L001_R2_001.fastq"
    ## [21] "F3D149_S215_L001_R1_001.fastq" "F3D149_S215_L001_R2_001.fastq"
    ## [23] "F3D150_S216_L001_R1_001.fastq" "F3D150_S216_L001_R2_001.fastq"
    ## [25] "F3D2_S190_L001_R1_001.fastq"   "F3D2_S190_L001_R2_001.fastq"  
    ## [27] "F3D3_S191_L001_R1_001.fastq"   "F3D3_S191_L001_R2_001.fastq"  
    ## [29] "F3D5_S193_L001_R1_001.fastq"   "F3D5_S193_L001_R2_001.fastq"  
    ## [31] "F3D6_S194_L001_R1_001.fastq"   "F3D6_S194_L001_R2_001.fastq"  
    ## [33] "F3D7_S195_L001_R1_001.fastq"   "F3D7_S195_L001_R2_001.fastq"  
    ## [35] "F3D8_S196_L001_R1_001.fastq"   "F3D8_S196_L001_R2_001.fastq"  
    ## [37] "F3D9_S197_L001_R1_001.fastq"   "F3D9_S197_L001_R2_001.fastq"  
    ## [39] "filtered"                      "HMP_MOCK.v35.fasta"           
    ## [41] "Mock_S280_L001_R1_001.fastq"   "Mock_S280_L001_R2_001.fastq"  
    ## [43] "mouse.dpw.metadata"            "mouse.time.design"            
    ## [45] "stability.batch"               "stability.files"

``` r
fnFs <- sort(list.files(miseq_path, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(miseq_path, pattern="_R2_001.fastq"))
sampleNames <- sapply(strsplit(fnFs, "_"), `[`, 1)
fnFs <- file.path(miseq_path, fnFs)
fnRs <- file.path(miseq_path, fnRs)
fnFs[1:3]
```

    ## [1] "/home/rstudio/MiSeq_SOP/F3D0_S188_L001_R1_001.fastq"  
    ## [2] "/home/rstudio/MiSeq_SOP/F3D1_S189_L001_R1_001.fastq"  
    ## [3] "/home/rstudio/MiSeq_SOP/F3D141_S207_L001_R1_001.fastq"

``` r
plotQualityProfile(fnFs[1:3])
```

    ## Warning: The `<scale>` argument of `guides()` cannot be `FALSE`. Use "none" instead as
    ## of ggplot2 3.3.4.
    ## ℹ The deprecated feature was likely used in the dada2 package.
    ##   Please report the issue at <]8;;https://github.com/benjjneb/dada2/issueshttps://github.com/benjjneb/dada2/issues]8;;>.

![](CC1-V2_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
plotQualityProfile(fnRs[1:3])
```

![](CC1-V2_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
filt_path <- file.path(miseq_path, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))
```

filtration qualité -\> trunclen après 240 sur read 1 et 160 sur read 2
(réduction longueur des seq comme devient mauvaise qualité -\> attention
laisser un chevauchement pour l’alignement) maxN=0 on ne veut pas de N
comme ça fait trop baisser la qualité / truncQ =2 en dessous de Q20 on
ne prend pas

``` r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)
```

    ##                               reads.in reads.out
    ## F3D0_S188_L001_R1_001.fastq       7793      7113
    ## F3D1_S189_L001_R1_001.fastq       5869      5299
    ## F3D141_S207_L001_R1_001.fastq     5958      5463
    ## F3D142_S208_L001_R1_001.fastq     3183      2914
    ## F3D143_S209_L001_R1_001.fastq     3178      2941
    ## F3D144_S210_L001_R1_001.fastq     4827      4312

Dereplication -\> ne pas garder en mémoire les séquences qui sont
plusieurs fois dans le fichier

``` r
derepFs <- derepFastq(filtFs, verbose=TRUE)
```

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D0_F_filt.fastq.gz

    ## Encountered 1979 unique sequences from 7113 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D1_F_filt.fastq.gz

    ## Encountered 1639 unique sequences from 5299 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D141_F_filt.fastq.gz

    ## Encountered 1477 unique sequences from 5463 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D142_F_filt.fastq.gz

    ## Encountered 904 unique sequences from 2914 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D143_F_filt.fastq.gz

    ## Encountered 939 unique sequences from 2941 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D144_F_filt.fastq.gz

    ## Encountered 1267 unique sequences from 4312 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D145_F_filt.fastq.gz

    ## Encountered 1756 unique sequences from 6741 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D146_F_filt.fastq.gz

    ## Encountered 1438 unique sequences from 4560 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D147_F_filt.fastq.gz

    ## Encountered 3590 unique sequences from 15637 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D148_F_filt.fastq.gz

    ## Encountered 2762 unique sequences from 11413 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D149_F_filt.fastq.gz

    ## Encountered 3021 unique sequences from 12017 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D150_F_filt.fastq.gz

    ## Encountered 1566 unique sequences from 5032 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D2_F_filt.fastq.gz

    ## Encountered 3707 unique sequences from 18075 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D3_F_filt.fastq.gz

    ## Encountered 1479 unique sequences from 6250 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D5_F_filt.fastq.gz

    ## Encountered 1195 unique sequences from 4052 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D6_F_filt.fastq.gz

    ## Encountered 1832 unique sequences from 7369 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D7_F_filt.fastq.gz

    ## Encountered 1183 unique sequences from 4765 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D8_F_filt.fastq.gz

    ## Encountered 1382 unique sequences from 4871 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D9_F_filt.fastq.gz

    ## Encountered 1709 unique sequences from 6504 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/Mock_F_filt.fastq.gz

    ## Encountered 897 unique sequences from 4314 total sequences read.

``` r
derepRs <- derepFastq(filtRs, verbose=TRUE)
```

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D0_R_filt.fastq.gz

    ## Encountered 1660 unique sequences from 7113 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D1_R_filt.fastq.gz

    ## Encountered 1349 unique sequences from 5299 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D141_R_filt.fastq.gz

    ## Encountered 1335 unique sequences from 5463 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D142_R_filt.fastq.gz

    ## Encountered 853 unique sequences from 2914 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D143_R_filt.fastq.gz

    ## Encountered 880 unique sequences from 2941 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D144_R_filt.fastq.gz

    ## Encountered 1286 unique sequences from 4312 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D145_R_filt.fastq.gz

    ## Encountered 1803 unique sequences from 6741 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D146_R_filt.fastq.gz

    ## Encountered 1265 unique sequences from 4560 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D147_R_filt.fastq.gz

    ## Encountered 3414 unique sequences from 15637 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D148_R_filt.fastq.gz

    ## Encountered 2522 unique sequences from 11413 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D149_R_filt.fastq.gz

    ## Encountered 2771 unique sequences from 12017 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D150_R_filt.fastq.gz

    ## Encountered 1415 unique sequences from 5032 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D2_R_filt.fastq.gz

    ## Encountered 3290 unique sequences from 18075 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D3_R_filt.fastq.gz

    ## Encountered 1390 unique sequences from 6250 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D5_R_filt.fastq.gz

    ## Encountered 1134 unique sequences from 4052 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D6_R_filt.fastq.gz

    ## Encountered 1635 unique sequences from 7369 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D7_R_filt.fastq.gz

    ## Encountered 1084 unique sequences from 4765 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D8_R_filt.fastq.gz

    ## Encountered 1161 unique sequences from 4871 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/F3D9_R_filt.fastq.gz

    ## Encountered 1502 unique sequences from 6504 total sequences read.

    ## Dereplicating sequence entries in Fastq file: /home/rstudio/MiSeq_SOP/filtered/Mock_R_filt.fastq.gz

    ## Encountered 732 unique sequences from 4314 total sequences read.

``` r
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames
```

création d’un modèle d’erreurs

``` r
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 33514080 total bases in 139642 reads from 20 samples will be used for learning the error rates.

``` r
errR <- learnErrors(filtRs, multithread=TRUE)
```

    ## 22342720 total bases in 139642 reads from 20 samples will be used for learning the error rates.

T2A = T donne A en bas, les q scores pour une position à score Q10 la
probabilité d’avoir un T donne A est :

``` r
plotErrors(errF)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](CC1-V2_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
plotErrors(errR)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis

![](CC1-V2_files/figure-gfm/unnamed-chunk-10-2.png)<!-- -->

correction des erreurs

``` r
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1979 unique sequences.
    ## Sample 2 - 5299 reads in 1639 unique sequences.
    ## Sample 3 - 5463 reads in 1477 unique sequences.
    ## Sample 4 - 2914 reads in 904 unique sequences.
    ## Sample 5 - 2941 reads in 939 unique sequences.
    ## Sample 6 - 4312 reads in 1267 unique sequences.
    ## Sample 7 - 6741 reads in 1756 unique sequences.
    ## Sample 8 - 4560 reads in 1438 unique sequences.
    ## Sample 9 - 15637 reads in 3590 unique sequences.
    ## Sample 10 - 11413 reads in 2762 unique sequences.
    ## Sample 11 - 12017 reads in 3021 unique sequences.
    ## Sample 12 - 5032 reads in 1566 unique sequences.
    ## Sample 13 - 18075 reads in 3707 unique sequences.
    ## Sample 14 - 6250 reads in 1479 unique sequences.
    ## Sample 15 - 4052 reads in 1195 unique sequences.
    ## Sample 16 - 7369 reads in 1832 unique sequences.
    ## Sample 17 - 4765 reads in 1183 unique sequences.
    ## Sample 18 - 4871 reads in 1382 unique sequences.
    ## Sample 19 - 6504 reads in 1709 unique sequences.
    ## Sample 20 - 4314 reads in 897 unique sequences.

``` r
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
```

    ## Sample 1 - 7113 reads in 1660 unique sequences.
    ## Sample 2 - 5299 reads in 1349 unique sequences.
    ## Sample 3 - 5463 reads in 1335 unique sequences.
    ## Sample 4 - 2914 reads in 853 unique sequences.
    ## Sample 5 - 2941 reads in 880 unique sequences.
    ## Sample 6 - 4312 reads in 1286 unique sequences.
    ## Sample 7 - 6741 reads in 1803 unique sequences.
    ## Sample 8 - 4560 reads in 1265 unique sequences.
    ## Sample 9 - 15637 reads in 3414 unique sequences.
    ## Sample 10 - 11413 reads in 2522 unique sequences.
    ## Sample 11 - 12017 reads in 2771 unique sequences.
    ## Sample 12 - 5032 reads in 1415 unique sequences.
    ## Sample 13 - 18075 reads in 3290 unique sequences.
    ## Sample 14 - 6250 reads in 1390 unique sequences.
    ## Sample 15 - 4052 reads in 1134 unique sequences.
    ## Sample 16 - 7369 reads in 1635 unique sequences.
    ## Sample 17 - 4765 reads in 1084 unique sequences.
    ## Sample 18 - 4871 reads in 1161 unique sequences.
    ## Sample 19 - 6504 reads in 1502 unique sequences.
    ## Sample 20 - 4314 reads in 732 unique sequences.

on prend les seq F et R dadaFs -\> corrigé derepFs -\> combien on avait
de seq (avant derepetition) ‘Mock’ commmunauté artificielle pour vérrif
quelque chose, on ne la veut pas

``` r
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
#Makesequencetable = table d'observations
seqtabAll <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])
table(nchar(getSequences(seqtabAll)))
```

    ## 
    ## 251 252 253 254 255 
    ##   1  85 186   5   2

un seq à 251 caractères etc..

enlever les chimères -\> ADNPol part puis cycle suivant fixation zone
conservée sur l’ADN d’une autre bact -\> résultat brin néosynthétisé mi
bact1 mi bact2

``` r
seqtabNoC <- removeBimeraDenovo(seqtabAll)
```

# Anotation taxonomique

on utilise RDP

``` bash
cd ~
wget https://zenodo.org/record/4310151/files/rdp_train_set_18.fa.gz
```

    ## --2022-11-09 09:45:12--  https://zenodo.org/record/4310151/files/rdp_train_set_18.fa.gz
    ## Resolving zenodo.org (zenodo.org)... 188.184.117.155
    ## Connecting to zenodo.org (zenodo.org)|188.184.117.155|:443... connected.
    ## HTTP request sent, awaiting response... 200 OK
    ## Length: 5760850 (5.5M) [application/octet-stream]
    ## Saving to: ‘rdp_train_set_18.fa.gz.3’
    ## 
    ##      0K .......... .......... .......... .......... ..........  0% 13.2M 0s
    ##     50K .......... .......... .......... .......... ..........  1% 15.5M 0s
    ##    100K .......... .......... .......... .......... ..........  2% 12.1M 0s
    ##    150K .......... .......... .......... .......... ..........  3% 12.6M 0s
    ##    200K .......... .......... .......... .......... ..........  4% 89.0M 0s
    ##    250K .......... .......... .......... .......... ..........  5% 14.0M 0s
    ##    300K .......... .......... .......... .......... ..........  6% 14.0M 0s
    ##    350K .......... .......... .......... .......... ..........  7% 17.8M 0s
    ##    400K .......... .......... .......... .......... ..........  7% 37.5M 0s
    ##    450K .......... .......... .......... .......... ..........  8% 16.2M 0s
    ##    500K .......... .......... .......... .......... ..........  9% 37.7M 0s
    ##    550K .......... .......... .......... .......... .......... 10% 10.0M 0s
    ##    600K .......... .......... .......... .......... .......... 11% 16.3M 0s
    ##    650K .......... .......... .......... .......... .......... 12% 20.2M 0s
    ##    700K .......... .......... .......... .......... .......... 13% 15.7M 0s
    ##    750K .......... .......... .......... .......... .......... 14% 21.0M 0s
    ##    800K .......... .......... .......... .......... .......... 15% 27.0M 0s
    ##    850K .......... .......... .......... .......... .......... 15% 26.5M 0s
    ##    900K .......... .......... .......... .......... .......... 16% 12.9M 0s
    ##    950K .......... .......... .......... .......... .......... 17% 20.6M 0s
    ##   1000K .......... .......... .......... .......... .......... 18% 15.3M 0s
    ##   1050K .......... .......... .......... .......... .......... 19% 29.1M 0s
    ##   1100K .......... .......... .......... .......... .......... 20% 8.33M 0s
    ##   1150K .......... .......... .......... .......... .......... 21% 77.3M 0s
    ##   1200K .......... .......... .......... .......... .......... 22% 12.3M 0s
    ##   1250K .......... .......... .......... .......... .......... 23% 41.5M 0s
    ##   1300K .......... .......... .......... .......... .......... 23% 13.4M 0s
    ##   1350K .......... .......... .......... .......... .......... 24% 97.7M 0s
    ##   1400K .......... .......... .......... .......... .......... 25% 15.8M 0s
    ##   1450K .......... .......... .......... .......... .......... 26% 12.3M 0s
    ##   1500K .......... .......... .......... .......... .......... 27%  147M 0s
    ##   1550K .......... .......... .......... .......... .......... 28% 15.5M 0s
    ##   1600K .......... .......... .......... .......... .......... 29%  121M 0s
    ##   1650K .......... .......... .......... .......... .......... 30% 14.5M 0s
    ##   1700K .......... .......... .......... .......... .......... 31% 62.4M 0s
    ##   1750K .......... .......... .......... .......... .......... 31% 14.5M 0s
    ##   1800K .......... .......... .......... .......... .......... 32% 16.5M 0s
    ##   1850K .......... .......... .......... .......... .......... 33% 71.9M 0s
    ##   1900K .......... .......... .......... .......... .......... 34% 16.9M 0s
    ##   1950K .......... .......... .......... .......... .......... 35% 46.0M 0s
    ##   2000K .......... .......... .......... .......... .......... 36% 19.8M 0s
    ##   2050K .......... .......... .......... .......... .......... 37% 41.9M 0s
    ##   2100K .......... .......... .......... .......... .......... 38% 17.9M 0s
    ##   2150K .......... .......... .......... .......... .......... 39% 42.0M 0s
    ##   2200K .......... .......... .......... .......... .......... 39% 19.4M 0s
    ##   2250K .......... .......... .......... .......... .......... 40% 37.3M 0s
    ##   2300K .......... .......... .......... .......... .......... 41% 17.8M 0s
    ##   2350K .......... .......... .......... .......... .......... 42% 11.2M 0s
    ##   2400K .......... .......... .......... .......... .......... 43% 68.3M 0s
    ##   2450K .......... .......... .......... .......... .......... 44% 19.3M 0s
    ##   2500K .......... .......... .......... .......... .......... 45% 40.0M 0s
    ##   2550K .......... .......... .......... .......... .......... 46% 15.1M 0s
    ##   2600K .......... .......... .......... .......... .......... 47% 16.7M 0s
    ##   2650K .......... .......... .......... .......... .......... 47% 87.5M 0s
    ##   2700K .......... .......... .......... .......... .......... 48% 15.4M 0s
    ##   2750K .......... .......... .......... .......... .......... 49% 13.2M 0s
    ##   2800K .......... .......... .......... .......... .......... 50% 77.6M 0s
    ##   2850K .......... .......... .......... .......... .......... 51% 14.7M 0s
    ##   2900K .......... .......... .......... .......... .......... 52% 53.7M 0s
    ##   2950K .......... .......... .......... .......... .......... 53% 11.9M 0s
    ##   3000K .......... .......... .......... .......... .......... 54% 34.3M 0s
    ##   3050K .......... .......... .......... .......... .......... 55% 9.51M 0s
    ##   3100K .......... .......... .......... .......... .......... 55% 54.8M 0s
    ##   3150K .......... .......... .......... .......... .......... 56% 9.55M 0s
    ##   3200K .......... .......... .......... .......... .......... 57% 65.4M 0s
    ##   3250K .......... .......... .......... .......... .......... 58% 14.3M 0s
    ##   3300K .......... .......... .......... .......... .......... 59% 61.9M 0s
    ##   3350K .......... .......... .......... .......... .......... 60% 15.2M 0s
    ##   3400K .......... .......... .......... .......... .......... 61% 74.2M 0s
    ##   3450K .......... .......... .......... .......... .......... 62% 15.5M 0s
    ##   3500K .......... .......... .......... .......... .......... 63% 75.2M 0s
    ##   3550K .......... .......... .......... .......... .......... 63% 17.9M 0s
    ##   3600K .......... .......... .......... .......... .......... 64% 14.7M 0s
    ##   3650K .......... .......... .......... .......... .......... 65% 90.7M 0s
    ##   3700K .......... .......... .......... .......... .......... 66% 18.3M 0s
    ##   3750K .......... .......... .......... .......... .......... 67% 16.8M 0s
    ##   3800K .......... .......... .......... .......... .......... 68% 15.5M 0s
    ##   3850K .......... .......... .......... .......... .......... 69% 63.2M 0s
    ##   3900K .......... .......... .......... .......... .......... 70% 13.8M 0s
    ##   3950K .......... .......... .......... .......... .......... 71% 56.2M 0s
    ##   4000K .......... .......... .......... .......... .......... 71% 16.1M 0s
    ##   4050K .......... .......... .......... .......... .......... 72% 89.4M 0s
    ##   4100K .......... .......... .......... .......... .......... 73% 16.3M 0s
    ##   4150K .......... .......... .......... .......... .......... 74% 84.0M 0s
    ##   4200K .......... .......... .......... .......... .......... 75% 18.1M 0s
    ##   4250K .......... .......... .......... .......... .......... 76% 50.5M 0s
    ##   4300K .......... .......... .......... .......... .......... 77% 15.6M 0s
    ##   4350K .......... .......... .......... .......... .......... 78% 60.6M 0s
    ##   4400K .......... .......... .......... .......... .......... 79% 17.9M 0s
    ##   4450K .......... .......... .......... .......... .......... 79% 76.5M 0s
    ##   4500K .......... .......... .......... .......... .......... 80% 16.1M 0s
    ##   4550K .......... .......... .......... .......... .......... 81% 82.2M 0s
    ##   4600K .......... .......... .......... .......... .......... 82% 16.1M 0s
    ##   4650K .......... .......... .......... .......... .......... 83% 16.5M 0s
    ##   4700K .......... .......... .......... .......... .......... 84% 49.8M 0s
    ##   4750K .......... .......... .......... .......... .......... 85% 19.2M 0s
    ##   4800K .......... .......... .......... .......... .......... 86% 50.8M 0s
    ##   4850K .......... .......... .......... .......... .......... 87% 21.0M 0s
    ##   4900K .......... .......... .......... .......... .......... 87% 17.2M 0s
    ##   4950K .......... .......... .......... .......... .......... 88% 57.7M 0s
    ##   5000K .......... .......... .......... .......... .......... 89% 12.2M 0s
    ##   5050K .......... .......... .......... .......... .......... 90%  127M 0s
    ##   5100K .......... .......... .......... .......... .......... 91% 14.5M 0s
    ##   5150K .......... .......... .......... .......... .......... 92% 14.0M 0s
    ##   5200K .......... .......... .......... .......... .......... 93% 78.3M 0s
    ##   5250K .......... .......... .......... .......... .......... 94% 13.6M 0s
    ##   5300K .......... .......... .......... .......... .......... 95% 96.5M 0s
    ##   5350K .......... .......... .......... .......... .......... 95% 12.6M 0s
    ##   5400K .......... .......... .......... .......... .......... 96% 34.4M 0s
    ##   5450K .......... .......... .......... .......... .......... 97% 17.1M 0s
    ##   5500K .......... .......... .......... .......... .......... 98% 67.1M 0s
    ##   5550K .......... .......... .......... .......... .......... 99% 7.32M 0s
    ##   5600K .......... .......... .....                           100%  126M=0.3s
    ## 
    ## 2022-11-09 09:45:13 (20.8 MB/s) - ‘rdp_train_set_18.fa.gz.3’ saved [5760850/5760850]

``` r
fastaRef <- "/home/rstudio/rdp_train_set_18.fa.gz"
taxTab <- assignTaxonomy(seqtabNoC, refFasta = fastaRef, multithread=TRUE)
unname(head(taxTab))
```

    ##      [,1]       [,2]            [,3]          [,4]            [,5]            
    ## [1,] "Bacteria" "Bacteroidetes" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [2,] "Bacteria" "Bacteroidetes" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [3,] "Bacteria" "Bacteroidetes" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [4,] "Bacteria" "Bacteroidetes" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ## [5,] "Bacteria" "Bacteroidetes" "Bacteroidia" "Bacteroidales" "Bacteroidaceae"
    ## [6,] "Bacteria" "Bacteroidetes" "Bacteroidia" "Bacteroidales" "Muribaculaceae"
    ##      [,6]         
    ## [1,] "Duncaniella"
    ## [2,] "Duncaniella"
    ## [3,] "Muribaculum"
    ## [4,] "Muribaculum"
    ## [5,] "Bacteroides"
    ## [6,] "Muribaculum"

``` r
seqs <- getSequences(seqtabNoC)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
```

``` r
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
        rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)
```

``` r
samdf <- read.csv("https://raw.githubusercontent.com/spholmes/F1000_workflow/master/data/MIMARKS_Data_combined.csv",header=TRUE)
samdf$SampleID <- paste0(gsub("00", "", samdf$host_subject_id), "D", samdf$age-21) #colonne ample id ds samdf
samdf <- samdf[!duplicated(samdf$SampleID),] # Remove dupicate entries for reverse reads
rownames(seqtabAll) <- gsub("124", "125", rownames(seqtabAll)) # Fix discrepancy #gsub = sub and replace -> manip de txt
all(rownames(seqtabAll) %in% samdf$SampleID) # TRUE
```

    ## [1] TRUE

``` r
samdf$SampleID
```

    ##   [1] "F3D0"   "F3D1"   "F3D2"   "F3D3"   "F3D4"   "F3D5"   "F3D6"   "F3D7"  
    ##   [9] "F3D8"   "F3D9"   "F3D11"  "F3D13"  "F3D15"  "F3D17"  "F3D19"  "F3D21" 
    ##  [17] "F3D25"  "F3D65"  "F3D125" "F3D141" "F3D142" "F3D143" "F3D144" "F3D145"
    ##  [25] "F3D146" "F3D147" "F3D148" "F3D149" "F3D150" "F3D364" "F4D0"   "F4D1"  
    ##  [33] "F4D2"   "F4D3"   "F4D4"   "F4D5"   "F4D6"   "F4D7"   "F4D8"   "F4D9"  
    ##  [41] "F4D11"  "F4D13"  "F4D15"  "F4D17"  "F4D19"  "F4D21"  "F4D25"  "F4D65" 
    ##  [49] "F4D125" "F4D141" "F4D142" "F4D143" "F4D144" "F4D145" "F4D146" "F4D147"
    ##  [57] "F4D148" "F4D149" "F4D150" "F5D0"   "F5D1"   "F5D2"   "F5D3"   "F5D4"  
    ##  [65] "F5D5"   "F5D6"   "F5D7"   "F5D8"   "F5D9"   "F5D11"  "F5D13"  "F5D15" 
    ##  [73] "F5D17"  "F5D19"  "F5D21"  "F5D25"  "F5D45"  "F5D65"  "F5D125" "F5D141"
    ##  [81] "F5D142" "F5D143" "F5D144" "F5D145" "F5D146" "F5D147" "F5D148" "F5D149"
    ##  [89] "F5D150" "F6D0"   "F6D1"   "F6D2"   "F6D3"   "F6D4"   "F6D5"   "F6D6"  
    ##  [97] "F6D7"   "F6D8"   "F6D9"   "F6D11"  "F6D13"  "F6D15"  "F6D17"  "F6D19" 
    ## [105] "F6D21"  "F6D25"  "F6D45"  "F6D65"  "F6D125" "F6D141" "F6D142" "F6D143"
    ## [113] "F6D144" "F6D145" "F6D146" "F6D147" "F6D148" "F6D149" "F6D150" "F7D0"  
    ## [121] "F7D1"   "F7D2"   "F7D3"   "F7D4"   "F7D5"   "F7D6"   "F7D7"   "F7D8"  
    ## [129] "F7D9"   "F7D11"  "F7D13"  "F7D15"  "F7D17"  "F7D19"  "F7D21"  "F7D25" 
    ## [137] "F7D45"  "F7D65"  "F7D125" "F7D141" "F7D142" "F7D143" "F7D144" "F7D145"
    ## [145] "F7D146" "F7D147" "F7D148" "F7D149" "F7D150" "F8D0"   "F8D1"   "F8D2"  
    ## [153] "F8D3"   "F8D4"   "F8D5"   "F8D6"   "F8D7"   "F8D8"   "F8D9"   "F8D11" 
    ## [161] "F8D13"  "F8D15"  "F8D17"  "F8D19"  "F8D21"  "F8D25"  "F8D45"  "F8D65" 
    ## [169] "F8D125" "F8D141" "F8D142" "F8D143" "F8D144" "F8D145" "F8D146" "F8D147"
    ## [177] "F8D148" "F8D149" "F8D150" "NADNA"  "M1D0"   "M1D1"   "M1D2"   "M1D3"  
    ## [185] "M1D4"   "M1D5"   "M1D6"   "M1D7"   "M1D8"   "M1D9"   "M1D11"  "M1D13" 
    ## [193] "M1D15"  "M1D17"  "M1D19"  "M1D21"  "M1D25"  "M1D65"  "M1D125" "M1D141"
    ## [201] "M1D142" "M1D143" "M1D144" "M1D145" "M1D146" "M1D147" "M1D148" "M1D149"
    ## [209] "M1D150" "M1D364" "M2D0"   "M2D1"   "M2D2"   "M2D3"   "M2D4"   "M2D5"  
    ## [217] "M2D6"   "M2D7"   "M2D8"   "M2D9"   "M2D11"  "M2D13"  "M2D15"  "M2D17" 
    ## [225] "M2D19"  "M2D21"  "M2D25"  "M2D65"  "M2D125" "M2D141" "M2D142" "M2D143"
    ## [233] "M2D144" "M2D145" "M2D146" "M2D147" "M2D148" "M2D149" "M2D150" "M2D364"
    ## [241] "M3D0"   "M3D1"   "M3D2"   "M3D3"   "M3D4"   "M3D5"   "M3D6"   "M3D7"  
    ## [249] "M3D8"   "M3D9"   "M3D11"  "M3D13"  "M3D15"  "M3D17"  "M3D19"  "M3D21" 
    ## [257] "M3D25"  "M3D45"  "M3D65"  "M3D125" "M3D141" "M3D142" "M3D143" "M3D144"
    ## [265] "M3D145" "M3D146" "M3D147" "M3D148" "M3D149" "M3D150" "M4D0"   "M4D1"  
    ## [273] "M4D2"   "M4D3"   "M4D4"   "M4D5"   "M4D6"   "M4D7"   "M4D8"   "M4D9"  
    ## [281] "M4D11"  "M4D13"  "M4D15"  "M4D17"  "M4D19"  "M4D21"  "M4D25"  "M4D45" 
    ## [289] "M4D65"  "M4D125" "M4D141" "M4D142" "M4D143" "M4D144" "M4D145" "M4D146"
    ## [297] "M4D147" "M4D148" "M4D149" "M4D150" "M5D0"   "M5D1"   "M5D2"   "M5D3"  
    ## [305] "M5D4"   "M5D5"   "M5D6"   "M5D7"   "M5D8"   "M5D9"   "M5D11"  "M5D13" 
    ## [313] "M5D15"  "M5D17"  "M5D19"  "M5D21"  "M5D25"  "M5D45"  "M5D65"  "M5D125"
    ## [321] "M5D141" "M5D142" "M5D143" "M5D144" "M5D145" "M5D146" "M5D147" "M5D148"
    ## [329] "M5D149" "M5D150" "M6D0"   "M6D1"   "M6D2"   "M6D3"   "M6D4"   "M6D5"  
    ## [337] "M6D6"   "M6D7"   "M6D8"   "M6D9"   "M6D11"  "M6D13"  "M6D15"  "M6D17" 
    ## [345] "M6D19"  "M6D21"  "M6D25"  "M6D45"  "M6D65"  "M6D125" "M6D141" "M6D142"
    ## [353] "M6D143" "M6D144" "M6D145" "M6D146" "M6D147" "M6D148" "M6D149" "M6D150"
    ## [361] "F6D155" "M4D155" "F5D165" "F6D165" "F7D175" "F8D175" "M3D175" "M4D175"
    ## [369] "M5D175" "M6D175" "F4D302" "F5D364" "F6D364" "F7D364" "F8D364" "M3D364"
    ## [377] "M4D364" "M5D364" "M6D364"

``` r
rownames(samdf) <- samdf$SampleID
keep.cols <- c("collection_date", "biome", "target_gene", "target_subfragment",
"host_common_name", "host_subject_id", "age", "sex", "body_product", "tot_mass",
"diet", "family_relationship", "genotype", "SampleID") 
samdf <- samdf[rownames(seqtabAll), keep.cols]

ps <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE), # on crée ps, un objet phyloseq -> objet phyloseq -> OTU table (chaque seq d'échantillon) -> accès ps@otu_table + Sample Data -> table avec echantillons, male femelle etc + TaxTable -> table de taxonomie avec seq phyla classe genre sp... + PhyThree -> arbre phylogen  ---> chacun avec sa structure diff
               sample_data(samdf), 
               tax_table(taxTab),phy_tree(fitGTR$tree))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 218 taxa and 19 samples ]
    ## sample_data() Sample Data:       [ 19 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 218 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 218 tips and 216 internal nodes ]

``` r
ps_connect <-url("https://raw.githubusercontent.com/spholmes/F1000_workflow/master/data/ps.rds") #ps - rds
ps = readRDS(ps_connect)
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 389 taxa and 360 samples ]
    ## sample_data() Sample Data:       [ 360 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 389 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 389 tips and 387 internal nodes ]

``` r
# Show available ranks in the dataset
rank_names(ps)
```

    ## [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"

``` r
test = ps@phy_tree
plot_tree(test)
```

![](CC1-V2_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
# Create table, number of features for each phyla
table(tax_table(ps)[, "Phylum"], exclude = NULL)
```

    ## 
    ##              Actinobacteria               Bacteroidetes 
    ##                          13                          23 
    ## Candidatus_Saccharibacteria   Cyanobacteria/Chloroplast 
    ##                           1                           4 
    ##         Deinococcus-Thermus                  Firmicutes 
    ##                           1                         327 
    ##                Fusobacteria              Proteobacteria 
    ##                           1                          11 
    ##                 Tenericutes             Verrucomicrobia 
    ##                           1                           1 
    ##                        <NA> 
    ##                           6

``` r
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized")) #on enleve NA ou uncharacterized ou vide

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))

plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
```

    ##                         Phylum         1     2
    ## 1               Actinobacteria 120.15385  1562
    ## 2                Bacteroidetes 265.52174  6107
    ## 3  Candidatus_Saccharibacteria 280.00000   280
    ## 4    Cyanobacteria/Chloroplast  64.25000   257
    ## 5          Deinococcus-Thermus  52.00000    52
    ## 6                   Firmicutes 179.24771 58614
    ## 7                 Fusobacteria   2.00000     2
    ## 8               Proteobacteria  59.09091   650
    ## 9                  Tenericutes 234.00000   234
    ## 10             Verrucomicrobia 104.00000   104

``` r
# Define phyla to filter
filterPhyla = c("Fusobacteria", "Deinococcus-Thermus")
# Filter entries with unidentified Phylum.
ps1 = subset_taxa(ps, !Phylum %in% filterPhyla)
ps1
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 381 taxa and 360 samples ]
    ## sample_data() Sample Data:       [ 360 samples by 14 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 381 taxa by 6 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 381 tips and 379 internal nodes ]

``` r
# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")
```

![](CC1-V2_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

``` r
# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps)
prevalenceThreshold
```

    ## [1] 18

``` r
# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps)
```

``` r
length(get_taxa_unique(ps2, taxonomic.rank = "Genus")) #Aglomerate taxa utilise la taxo pour regroup des seq qui ont le même rang taxonomique -> ts les seq ds le même Genre se mettent ds un groupe du Genre -> nouvelle unité écologique, taxa pas OTU
```

    ## [1] 49

``` r
ps3 = tax_glom(ps2, "Genus", NArm = TRUE)

h1 = 0.4
ps4 = tip_glom(ps2, h = h1)

multiPlotTitleTextSize = 15
p2tree = plot_tree(ps2, method = "treeonly",
                   ladderize = "left",
                   title = "Before Agglomeration") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p3tree = plot_tree(ps3, method = "treeonly",
                   ladderize = "left", title = "By Genus") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p4tree = plot_tree(ps4, method = "treeonly",
                   ladderize = "left", title = "By Height") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))

# group plots together
grid.arrange(nrow = 1, p2tree, p3tree, p4tree)
```

![](CC1-V2_files/figure-gfm/unnamed-chunk-27-1.png)<!-- --> 1 -\>
données brutes  
2 -\> regroupées par Genres 3 -\> regroupées par distance ds le premier
arbre

``` r
#Normalise le jeu de données -> ici normalisation de variances
plot_abundance = function(physeq,title = "",
                          Facet = "Order", Color = "Phylum"){
  # Arbitrary subset, based on Phylum, for plotting
  p1f = subset_taxa(physeq, Phylum %in% c("Firmicutes"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "sex",y = "Abundance",
                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}

# Transform to relative abundance. Save as new object.
ps3ra = transform_sample_counts(ps3, function(x){x / sum(x)}) # pour chaque cell on divise par la somme des cell

plotBefore = plot_abundance(ps3,"")
```

    ## Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
    ## ℹ Please use tidy evaluation ideoms with `aes()`

``` r
plotAfter = plot_abundance(ps3ra,"")
# Combine each plot into one graphic.
grid.arrange(nrow = 2,  plotBefore, plotAfter)
```

![](CC1-V2_files/figure-gfm/unnamed-chunk-28-1.png)<!-- --> avant et
après normalisation

``` r
psOrd = subset_taxa(ps3ra, Order == "Lactobacillales")
plot_abundance(psOrd, Facet = "Genus", Color = NULL)
```

![](CC1-V2_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

``` r
qplot(sample_data(ps)$age, geom = "histogram",binwidth=20) + xlab("age")
```

    ## Warning: `qplot()` was deprecated in ggplot2 3.4.0.

![](CC1-V2_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

``` r
qplot(log10(rowSums(otu_table(ps))),binwidth=0.2) +
  xlab("Logged counts-per-sample")
```

![](CC1-V2_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

``` r
sample_data(ps)$age_binned <- cut(sample_data(ps)$age,
                          breaks = c(0, 100, 200, 400))
levels(sample_data(ps)$age_binned) <- list(Young100="(0,100]", Mid100to200="(100,200]", Old200="(200,400]")
sample_data(ps)$family_relationship=gsub(" ","",sample_data(ps)$family_relationship)
pslog <- transform_sample_counts(ps, function(x) log(1 + x))
out.wuf.log <- ordinate(pslog, method = "MDS", distance = "wunifrac")
```

    ## Warning in UniFrac(physeq, weighted = TRUE, ...): Randomly assigning root as --
    ## GCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGACGGAATAGCAAGTCTGATGTGAAAACCCGGGGCTCAACTCCGGGAGTGCATTGGAAACTGTTGATCTAGAGTGTCGGAGAGGCAAGCGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAACACCAGTGGCGAAGGCGGCTTGCTGGACGATGACTGACGTTGAGGCTCGAAAGCGTGGGGAG
    ## -- in the phylogenetic tree in the data you provided.

``` r
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "age_binned") +
  labs(col = "Binned Age") +
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](CC1-V2_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

``` r
rel_abund <- t(apply(otu_table(ps), 1, function(x) x / sum(x)))
qplot(rel_abund[, 12], geom = "histogram",binwidth=0.05) +
  xlab("Relative abundance")
```

![](CC1-V2_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

``` r
outliers <- c("F5D165", "F6D165", "M3D175", "M4D175", "M5D175", "M6D175")
ps <- prune_samples(!(sample_names(ps) %in% outliers), ps)

which(!rowSums(otu_table(ps)) > 1000)
```

    ## F5D145 M1D149   M1D9 M2D125  M2D19 M3D148 M3D149   M3D3   M3D5   M3D8 
    ##     69    185    200    204    218    243    244    252    256    260

``` r
ps <- prune_samples(rowSums(otu_table(ps)) > 1000, ps)
pslog <- transform_sample_counts(ps, function(x) log(1 + x))

out.pcoa.log <- ordinate(pslog,  method = "MDS", distance = "bray")
evals <- out.pcoa.log$values[,1]
plot_ordination(pslog, out.pcoa.log, color = "age_binned",
                  shape = "family_relationship") +
  labs(col = "Binned Age", shape = "Litter")+
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](CC1-V2_files/figure-gfm/unnamed-chunk-34-1.png)<!-- -->

``` r
out.dpcoa.log <- ordinate(pslog, method = "DPCoA")
evals <- out.dpcoa.log$eig
plot_ordination(pslog, out.dpcoa.log, color = "age_binned", label= "SampleID",
                  shape = "family_relationship") +
  labs(col = "Binned Age", shape = "Litter")+
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](CC1-V2_files/figure-gfm/unnamed-chunk-35-1.png)<!-- -->

``` r
plot_ordination(pslog, out.dpcoa.log, type = "species", color = "Phylum") +
  coord_fixed(sqrt(evals[2] / evals[1]))
```

![](CC1-V2_files/figure-gfm/unnamed-chunk-36-1.png)<!-- -->

``` r
out.wuf.log <- ordinate(pslog, method = "PCoA", distance ="wunifrac")
```

    ## Warning in UniFrac(physeq, weighted = TRUE, ...): Randomly assigning root as --
    ## GCAAGCGTTGTCCGGAATTACTGGGTGTAAAGGGAGCGTAGGCGGGAAAGCAAGTTGAATGTTTAAGGTACGGGCTCAACCTGTACAAGCGTTCAAAACTGTTTTTCTTGAGTGGAGTAGAGGTAAGCGGAATTCCTAGTGTAGCGGTGAAATGCGTAAATATTAGGAGGAACACCAGTGGCGAAGGCGGCTTACTGGGCTTTAACTGACGCTGAGGCTCGAAAGCGTGGGTAGC
    ## -- in the phylogenetic tree in the data you provided.

``` r
evals <- out.wuf.log$values$Eigenvalues
plot_ordination(pslog, out.wuf.log, color = "age_binned",
                  shape = "family_relationship") +
  coord_fixed(sqrt(evals[2] / evals[1])) +
  labs(col = "Binned Age", shape = "Litter")
```

![](CC1-V2_files/figure-gfm/unnamed-chunk-37-1.png)<!-- -->

``` r
abund <- otu_table(pslog)
abund_ranks <- t(apply(abund, 1, rank))
abund_ranks <- abund_ranks - 329
abund_ranks[abund_ranks < 1] <- 1

library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following objects are masked from 'package:Biostrings':
    ## 
    ##     collapse, intersect, setdiff, setequal, union

    ## The following object is masked from 'package:GenomeInfoDb':
    ## 
    ##     intersect

    ## The following object is masked from 'package:XVector':
    ## 
    ##     slice

    ## The following objects are masked from 'package:IRanges':
    ## 
    ##     collapse, desc, intersect, setdiff, slice, union

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, intersect, rename, setdiff, setequal, union

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(reshape2)
abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
```

    ## Joining, by = c("Var1", "Var2")

``` r
colnames(abund_df) <- c("sample", "seq", "abund", "rank")

abund_df <- melt(abund, value.name = "abund") %>%
  left_join(melt(abund_ranks, value.name = "rank"))
```

    ## Joining, by = c("Var1", "Var2")

``` r
colnames(abund_df) <- c("sample", "seq", "abund", "rank")

sample_ix <- sample(1:nrow(abund_df), 8)
ggplot(abund_df %>%
         filter(sample %in% abund_df$sample[sample_ix])) +
  geom_point(aes(x = abund, y = rank, col = sample),
             position = position_jitter(width = 0.2), size = 1.5) +
  labs(x = "Abundance", y = "Thresholded rank") +
  scale_color_brewer(palette = "Set2")
```

![](CC1-V2_files/figure-gfm/unnamed-chunk-38-1.png)<!-- -->

``` r
library(ade4)
```

    ## 
    ## Attaching package: 'ade4'

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     score

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     score

``` r
ranks_pca <- dudi.pca(abund_ranks, scannf = F, nf = 3)
row_scores <- data.frame(li = ranks_pca$li,
                         SampleID = rownames(abund_ranks))
col_scores <- data.frame(co = ranks_pca$co,
                         seq = colnames(abund_ranks))
tax <- tax_table(ps) %>%
  data.frame(stringsAsFactors = FALSE)
tax$seq <- rownames(tax)
main_orders <- c("Clostridiales", "Bacteroidales", "Lactobacillales",
                 "Coriobacteriales")
tax$Order[!(tax$Order %in% main_orders)] <- "Other"
tax$Order <- factor(tax$Order, levels = c(main_orders, "Other"))
tax$otu_id <- seq_len(ncol(otu_table(ps)))
row_scores <- row_scores %>%
  left_join(sample_data(pslog))
```

    ## Joining, by = "SampleID"

``` r
col_scores <- col_scores %>%
  left_join(tax)
```

    ## Joining, by = "seq"

``` r
evals_prop <- 100 * (ranks_pca$eig / sum(ranks_pca$eig))
ggplot() +
  geom_point(data = row_scores, aes(x = li.Axis1, y = li.Axis2), shape = 2) +
  geom_point(data = col_scores, aes(x = 25 * co.Comp1, y = 25 * co.Comp2, col = Order),
             size = .3, alpha = 0.6) +
  scale_color_brewer(palette = "Set2") +
  facet_grid(~ age_binned) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
       y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
  coord_fixed(sqrt(ranks_pca$eig[2] / ranks_pca$eig[1])) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
```

![](CC1-V2_files/figure-gfm/unnamed-chunk-39-1.png)<!-- -->

``` r
ps_ccpna <- ordinate(pslog, "CCA", formula = pslog ~ age_binned + family_relationship)

library(ggrepel)
ps_scores <- vegan::scores(ps_ccpna)
sites <- data.frame(ps_scores$sites)
sites$SampleID <- rownames(sites)
sites <- sites %>%
  left_join(sample_data(ps))
```

    ## Joining, by = "SampleID"

``` r
species <- data.frame(ps_scores$species)
species$otu_id <- seq_along(colnames(otu_table(ps)))
species <- species %>%
  left_join(tax)
```

    ## Joining, by = "otu_id"

``` r
evals_prop <- 100 * ps_ccpna$CCA$eig[1:2] / sum(ps_ccpna$CA$eig)
ggplot() +
  geom_point(data = sites, aes(x = CCA1, y = CCA2), shape = 2, alpha = 0.5) +
  geom_point(data = species, aes(x = CCA1, y = CCA2, col = Order), size = 0.5) +
  geom_text_repel(data = species %>% filter(CCA2 < -2),
                    aes(x = CCA1, y = CCA2, label = otu_id),
            size = 1.5, segment.size = 0.1) +
  facet_grid(. ~ family_relationship) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
        y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
  scale_color_brewer(palette = "Set2") +
  coord_fixed(sqrt(ps_ccpna$CCA$eig[2] / ps_ccpna$CCA$eig[1])*0.45   ) +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))
```

    ## Warning: ggrepel: 9 unlabeled data points (too many overlaps). Consider increasing max.overlaps
    ## ggrepel: 9 unlabeled data points (too many overlaps). Consider increasing max.overlaps

![](CC1-V2_files/figure-gfm/unnamed-chunk-40-1.png)<!-- -->
