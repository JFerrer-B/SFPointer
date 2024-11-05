TUTORIAL. SFPointer: A Systematic Identification of RBPs Driving Aberrant Splicing in Cancer
================
César Lobato, Marian Gimeno, Ángel Rubio and Juan A. Ferrer-Bonsoms

*Tecnun (University of Navarra), Paseo Manuel Lardizábal 15, 20018 San Sebastián, SPAIN*

<br />

This is the tutorial for the publication: *"SFPointer: A Systematic Identification of RBPs Driving Aberrant Splicing in Cancer"*

In this tutorial, we present the methodology to predict RBPs based on the splicing events predicted by EventPointer. For this we use a combination of CLIP experiments in POSTAR3 with transcriptomic data (See the main manuscript for more details).

Here, we show the pipeline using a experiment in which TDP43 specific splicing factor is knocked-down.

The SFPointer methodology has different stages:
- Apply EventPointer (EP) to find differentially spliced ​​events between conditions.
- Apply SFPointer to find enriched RBPs.


Apply EP
---------------

For more details about EP go to https://www.bioconductor.org/packages/release/bioc/vignettes/EventPointer/inst/doc/EventPointer.html.

To apply EP we need the following input files:
- Reference transcriptome: in this case we used Gencode 24.
- Pseudoalignnet output. In this case, we used Kallisto.

As described in the EP vignette, there are several processes:
- Find the splicing events (this has to be done only once for each reference transcriptome)
- Compute the PSI.
- Apply the bootstrap method to find differentially spliced ​​events.

The following code chuck shows how to find the splicing events using the function *EventDetection_transcriptome*:

``` r
library(parallel)
library(EventPointer)
# only need to be run once per reference transcriptome
Eventsxisof_gc24 <- EventDetection_transcriptome(inputFile = "path to transcriptome/gencode24.gtf",
                                                 Transcriptome = "grch38",
                                                 Pathtxt = "./EP_events/",
                                                 cores = 16)
```

Given the output of *EventDetection_transcriptome* and the results obtained from the pseudoalignment, we can compute the PSI values:

``` r
pathfiles <- dir("path to kallisto output",full.names = TRUE)
RNA_seq <- getbootstrapdata(PathSamples = pathfiles,type = "kallisto")

PSI_evexp <- GetPSI_FromTranRef(Samples = RNA_seq,
                                PathsxTranscript = Eventsxisof_gc24,
                                Bootstrap = T,
                                Filter = T)
```

Finally, we can compute the boostrap statistic test:

``` r
Dmatrix <- cbind(c(1,1,1,0,0,0),
                 c(0,0,0,1,1,1))

Cmatrix <- as.matrix(limma::makeContrasts(TDP43-KOTDP43, levels=c("TDP43", "KOTDP43")))
ResultBootstrap_kallisto_TDP43 <- EventPointer_Bootstraps(PSI = PSI_evexp$PSI,
                                                 Design = Dmatrix,
                                                 Contrast = Cmatrix,
                                                 cores = 6,
                                                 ram = 20,
                                                 nBootstraps = 1000,
                                                 UsePseudoAligBootstrap = 1,
                                                 Threshold = 0)
```


SF Pointer
------------------

To apply SFPointer we need the following input data:

- The result of Event Pointer.
- The ExS matrix.

We can load this information from the *data* directory.
``` r
load("./data/ResultBootstrap_kallisto_TDP43.RData")
load("./data/ExS.RData")
```

Now, we can load the corresponding functions of SFPointer:
``` r
source("./SF_pointer_functions/SF_Prediction_functions.R")
source("./SF_pointer_functions/aux_functions_SF_Prediction.r")
```

SFPointer offers 4 different methods for enrichment: 
- Fisher's exact test, 
- Poisson Binomial, 
- Wilcoxon test
- GSEA.

The first two need to set which events are differentialy spliced. For this, you can select the top n events with the lowest p-value (default n = 1000) or the events with the p.value lower than th (by default, th = 0.001). Furthermore, you can set the th using the FDR information returned by EventPointer.

``` r
maximumFDR = 0.1
events_index = which(ResultBootstrap_kallisto_TDP43$LocalFDR$qvalues <= maximumFDR)
pval_th = max(ResultBootstrap_kallisto_TDP43$LocalFDR$pvalues[events_index])
```

This is, setting a threshold of 0.0016 to the p-value to consider a event to be differentialy spliced, corresponds to a FDR of 0.1.

Now, we can apply the Fishers and the Poisson Binomial approach:

``` r
SF_tdp43_Fisher <- SF_Prediction(P_value_PSI = 
          ResultBootstrap_kallisto_TDP43$Pvalues,
          ExS = ExS, method = "Fisher",
          significance = pval_th)
head(SF_tdp43_Fisher[[1]])
```
         RBP nHits Pvalue_hyp_PSI     N     d    n   x qhyp_0.5       Fc
TDP43   TDP43  4287   5.354573e-50 67679  4287 1089 213       69 3.086957
EIF4G2 EIF4G2   957   1.221520e-26 67679   957 1089  71       15 4.733333
YTHDC2 YTHDC2  3218   8.798743e-26 67679  3218 1089 138       52 2.653846
DKC1     DKC1   910   1.264320e-23 67679   910 1089  65       14 4.642857
ALKBH1 ALKBH1  4676   1.403874e-22 67679  4676 1089 167       75 2.226667
FIP1L1 FIP1L1 36883   2.134727e-22 67679 36883 1089 748      593 1.261383

The table shows the top 6 ranked RBPS. As expected, TDP 43 is ranked as the most likely RBP to be driving the differences in the splicing patterns between the two conditions.

``` r
SF_tdp43_PoiBin <- SF_Prediction(P_value_PSI = 
          ResultBootstrap_kallisto_TDP43$Pvalues,
          ExS = ExS, method = "PoiBin",
          significance = pval_th)
head(SF_tdp43_PoiBin[[1]])
```
          RBP nHits Pvalue_hyp_PSI     N    d    n   x  qhyp_0.5       Fc
EIF4G2 EIF4G2   957   0.000000e+00 67679  957 1089  71  26.69602 2.659572
TDP43   TDP43  4287   0.000000e+00 67679 4287 1089 213 110.83233 1.921822
DKC1     DKC1   910   1.776357e-15 67679  910 1089  65  25.41671 2.557372
YTHDC2 YTHDC2  3218   3.452121e-10 67679 3218 1089 138  85.04442 1.622681
UNR       UNR  2149   1.337532e-08 67679 2149 1089  99  58.18057 1.701599
CDC40   CDC40  1827   1.562860e-07 67679 1827 1089  85  49.84773 1.705193

The table shows the top 6 ranked RBPS. As expected, TDP 43 is ranked with EIF4G2 as the most likely RBPs to be driving the differences in the splicing patterns between the two conditions.



SFPointer includes the Wilcoxon and the GSEA approach, which are non-parametric enrichment analysis. Although you do not need to set a th to perform enrichment, SFPointer gives you the possibility to only work with events that are significant. For these cases, we will consider significant from a p-value of 0.05.

``` r
SF_tdp43_Wilcoxon <- SF_Prediction(P_value_PSI = 
          ResultBootstrap_kallisto_TDP43$Pvalues,
          ExS = ExS, method = "Wilcoxon",
          significance = rep(0.05,4))
head(SF_tdp43_Wilcoxon[[1]]) 
```
          RBP   z_score Pvalue_hyp_PSI
TDP43   TDP43 -9.174647   2.265278e-20
POU5F1 POU5F1 -6.256238   1.971879e-10
RBFOX1 RBFOX1 -5.872722   2.143480e-09
ALKBH1 ALKBH1 -5.800011   3.315522e-09
PABPC1 PABPC1 -5.747185   4.537063e-09
CPSF4   CPSF4 -5.717762   5.396813e-09



``` r
SF_tdp43_Gsea <- SF_Prediction(P_value_PSI = 
          ResultBootstrap_kallisto_TDP43$Pvalues,
          ExS = ExS, method = "Gsea",
          significance = rep(0.05,4))
head(SF_tdp43_Gsea[[1]]) 
```
       RBP         pval         padj   log2err        ES      NES size
204  TDP43 4.936410e-21 1.199548e-18 1.1866510 0.2338165 5.442541  497
61  EIF4G2 4.692082e-11 5.700879e-09 0.8513391 0.2932383 3.949430  143
231 YTHDC2 2.334920e-10 1.891285e-08 0.8140358 0.1870098 3.783252  350
52    DKC1 6.983255e-10 3.686742e-08 0.8012156 0.2867855 3.734236  131
9   ALKBH1 7.585889e-10 3.686742e-08 0.8012156 0.1649492 3.771187  473
20   CDC40 2.750174e-08 1.113820e-06 0.7337620 0.2228016 3.444384  193

