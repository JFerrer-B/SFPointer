TUTORIAL. SFPointer: A Systematic Identification of RBPs Driving Aberrant Splicing in Cancer
================
César Lobato, Marian Gimeno, Ángel Rubio and Juan A. Ferrer-Bonsoms

*Tecnun (University of Navarra), Paseo Manuel Lardizábal 15, 20018 San Sebastián, SPAIN*

<br />

This is the tutorial for the publication: *"SFPointer: A Systematic Identification of RBPs Driving Aberrant Splicing in Cancer"*

In this tutorial, we present the methodology to predict RBPs based on the splicing events predicted by EventPointer. For this we use a combination of CLIP experiments in POSTAR3 with transcriptomic data (See the main manuscript for more details).

Here, we show the pipeline using a experiment in which TDP43 specific splicing factor is knocked-down.


Library and source
------------------



``` r
library(parallel)
library(EventPointer)
load("./data/ResultsBootsTrapTDP43.RData")
load("./ExS_Postar3.RData")
```
EventPointer event detection
------------------
``` r

######1. detect all possible events in reference transcriptome ######
Eventsxisof_grch38 <- EventDetection_transcriptome(inputFile = "./Homo_sapiens.GRCh38.gtf",Transcriptome = "grch38",
                                                       Pathtxt = "./EP_events/",cores = 16)

#load("../Paper_SF_pointer_2_0/info_gencode24/EventsXtrans_gc24.RData")
######2. get psi of all the events ######


pathfiles <- dir("./GSE136366/kallisto_output/",full.names = TRUE)
RNA_seq <- getbootstrapdata(PathSamples = pathfiles,type = "kallisto")
save(RNA_seq,file="RNA_seq_GSE136366.RData")

PSI_evexp <- GetPSI_FromTranRef(Samples = RNA_seq,
                                PathsxTranscript = Eventsxisof_grch38,
                                Bootstrap = T,
                                Filter = T,
                                Qn = 0.001)

######3. run statistic (both the ST and ML version) ######
Cmatrix <- as.matrix(limma::makeContrasts(TDP43-KOTDP43, levels=c("TDP43", "KOTDP43")))


Dmatrix <- cbind(c(1,1,1,0,0,0),
                 c(0,0,0,1,1,1))

ResultsBootsTrapTDP43 <- EventPointer_Bootstraps(PSI = PSI_evexp$PSI,
                                                 Design = Dmatrix,
                                                 Contrast = Cmatrix,
                                                 cores = 6,
                                                 ram = 20,
                                                 nBootstraps = 1000,
                                                 UsePseudoAligBootstrap = 1,
                                                 Threshold = 0)
```

SFPrediction
---------------

``` r
######1. Fisher test ######
hyper_pred_TDP43 <- EventPointer:::SF_Prediction(ResultsBootsTrapTDP43$Pvalues, ExS_Postar3, method = "Fisher")
which(hyper_pred_TDP43[[1]][order(hyper_pred_TDP43[[1]]$Pvalue_hyp_PSI,decreasing = F),]$RBP == "TDP43" )
save(hyper_pred_TDP43,file = "hyper_pred_TDP43.RData")

######1. Poisson Binomial test ######
poisonBi_pred_TDP43 <- EventPointer:::SF_Prediction(ResultsBootsTrapTDP43$Pvalues, ExS_Postar3, method = "PoiBin")
which(poisonBi_pred_TDP43[[1]][order(poisonBi_pred_TDP43[[1]]$Pvalue_hyp_PSI,decreasing = F),]$RBP == "TDP43" )
save(poisonBi_pred_TDP43,file = "poisonBi_pred_TDP43.RData")

######1. Wilcoxon test ######
Wilcoxon_pred_TDP43 <- EventPointer:::SF_Prediction(ResultsBootsTrapTDP43$Pvalues, ExS_Postar3,method = "Wilcoxon")
which(Wilcoxon_pred_TDP43[[1]][order(Wilcoxon_pred_TDP43[[1]]$z.score,decreasing = F),]$RBP == "TDP43" )
save(Wilcoxon_pred_TDP43,file = "Wilcoxon_pred_TDP43.RData")

######1. GSEA test ######
Gsea_pred_TDP43 <- EventPointer:::SF_Prediction(ResultsBootsTrapTDP43$Pvalues, ExS_Postar3, method = "Gsea")
which(Gsea_pred_TDP43[[1]][order(Gsea_pred_TDP43[[1]]$pval,decreasing = F),]$pathway == "TDP43" )
save(Gsea_pred_TDP43,file = "Gsea_pred_TDP43.RData")
save(ResultsBootsTrapTDP43, file="ResultsBootsTrapTDP43.RData")
```

![](Tutorial-gitlab_files/figure-markdown_github/unnamed-chunk-10-1.png)![](Tutorial-gitlab_files/figure-markdown_github/unnamed-chunk-10-2.png)![](Tutorial-gitlab_files/figure-markdown_github/unnamed-chunk-10-3.png)![](Tutorial-gitlab_files/figure-markdown_github/unnamed-chunk-10-4.png)![](Tutorial-gitlab_files/figure-markdown_github/unnamed-chunk-10-5.png)![](Tutorial-gitlab_files/figure-markdown_github/unnamed-chunk-10-6.png)![](Tutorial-gitlab_files/figure-markdown_github/unnamed-chunk-10-7.png)![](Tutorial-gitlab_files/figure-markdown_github/unnamed-chunk-10-8.png)![](Tutorial-gitlab_files/figure-markdown_github/unnamed-chunk-10-9.png)![](Tutorial-gitlab_files/figure-markdown_github/unnamed-chunk-10-10.png)![](Tutorial-gitlab_files/figure-markdown_github/unnamed-chunk-10-11.png)![](Tutorial-gitlab_files/figure-markdown_github/unnamed-chunk-10-12.png)
