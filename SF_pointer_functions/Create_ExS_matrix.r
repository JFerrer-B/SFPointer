
CreateExSmatrix <- function(POSTAR,
                            EventsRegions){
  
  GRseq3 <- EventsRegions
  
  # Load PSTAR3 peaks #####
  POSTAR_L <- split.data.frame(x = POSTAR, f = POSTAR$RBP)
  
  
  mySF <- unique(names(POSTAR_L))
  # mySF <- mySF[!mySF %in% c("AGO2MNASE", "HURMNASE", "PAN-ELAVL","RBP_occupancy")]
  nSF <- length(mySF)
  
  EventID <- unique(EventsRegions$EventID)
  
  ExS <- matrix(0, nrow = length(EventID), ncol = nSF)
  rownames(ExS) <- EventID
  colnames(ExS) <- mySF
  
  for(i in seq_len(nSF)){
    # i <- 1
    # i <- 5
    SF<- mySF[i]
    # cat(i, "/", nSF, " ",as.character(SF), "\n")
    jjx <- match(SF, names(POSTAR_L))
    if(!is.na(jjx)){
      # cat("holaaaa \n\n")
      peaks <- POSTAR_L[[jjx]]
      iD <- which(!(peaks$seqname %in% c(paste0("chr", c(1:22)), "chrX", "chrY")))
      if(length(iD) > 0){
        peaks <- peaks[-iD, ]
      }
      peaks_GR <- GRanges(peaks)
      seqlevels(peaks_GR) <- sub("chr", "", seqlevels(peaks_GR))
      peaks_GR<-reduce(peaks_GR)
      Overlaps <- findOverlaps(peaks_GR, GRseq3)
      EvMatch <- as.character(elementMetadata(GRseq3)$EventID[subjectHits(Overlaps)])
      
      if(length(EvMatch)>0){
        ExS[EvMatch,i] <- 1
      }
    }
  }
  
  return(ExS)    
  
}




Compute_Events_Regions <- function(EventsFound,
                                   SG_List,
                                   cores,
                                   nt_intron=300,
                                   nt_exon=100,
                                   remove_complex=TRUE){
  
  EventsFound$EventID <- paste0(EventsFound$GeneName,"_",EventsFound$EventNumber)
  #type of events
  typeA <- c("Cassette Exon",
             "Retained Intron",
             "Alternative 3' Splice Site",
             "Alternative 5' Splice Site")
  
  
  # A cassette Event, for example, has de following structure.
  # The key points in this case are the points 1,2,3,4
  # REF5' |||||-------------||||||||------------||||| REF3'
  # REF5' |||||---------------------------------||||| REF3'
  #          a   b        c  d     e  f        g   h
  
  # retained intron:
  # REF5' ||||||||||||||||||||||||||||||||||||||||||||||| REF3'
  # REF5' |||||||||||||--------------------|||||||||||||| REF3'
  #                  c  d                e   f
  
  # alternative 3':
  # REF5' |||||-------------------------|||||||| REF3'
  # REF5' |||||---------------|||||||||||||||||| REF3'
  #          a  b          c   d       e  f
  
  # alternative 5':
  # REF5' ||||||||||||||||-------------||||||||| REF3'
  # REF5' |||||------------------------||||||||| REF3'
  #          c  d       e   f       g   h
  
  # alternative first:
  # REF5'                    ||||||||-----------||||| REF3'
  # REF5' |||||---------------------------------||||| REF3'
  #                       c  d     e  f        g   h
  
  # alternative Last:
  # REF5' |||||-----------|||||||||| REF3'
  # REF5' |||||---------------------------------||||| REF3'
  #         a   b        c  d     e  f        
  
  # compute Events_Regions ####
  
  GRseq3 <- callGRseq_parallel(EventsFound = EventsFound,
                               SG_List = SG_List,
                               cores = cores,
                               typeA = typeA,
                               nt_intron=nt_intron,
                               nt_exon=nt_exon,
                               remove_complex)
  
  return(GRseq3)
  
}




callGRseq_parallel <- function(EventsFound,SG_List,cores,typeA,nt_intron,nt_exon,remove_complex){
  
  
  #check if reclasification of events exists:
  # if(!is.null(EventsFound$EventType_new)){
  #   wwx <- which(EventsFound$EventType_new == "")
  #   if(length(wwx)>0) EventsFound$EventType_new[wwx] <- EventsFound$EventType[wwx]
  #   EventsFound$EventType <- EventsFound$EventType_new
  # }
  
  if(remove_complex){
    cat("\n removing complex events. By default, complex events are not consider for 
        this analysis:\n")
    rrx <- which(EventsFound$EventType == "Complex Event")
    EventsFound <- EventsFound[-rrx,,drop=FALSE]
    cat("removing ",length(rrx)," complex events from analysis.")
  }
  # table(EventsFound$EventType_new)
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  # unique(EventsFound$EventType)
  # which(EventsFound$EventType=="Cassette Exon")[1]
  # which(EventsFound$EventType=="Alternative 3' Splice Site")[1]
  # son5prima <- which(EventsFound$EventType=="Alternative 5' Splice Site")
  # mygene <- EventsFound$GeneName[son5prima]
  # for (m in 1:length(mygene)) {
  #   jjx <- which(names(SG_List)==mygene[m])
  #   myEdges <- SG_List[[jjx]]$Edges
  #   if(myEdges$Strand[1]=="-") break
  # }
  # son5prima[m]
  # 
  # 
  # which(EventsFound$EventType=="Mutually Exclusive Exons")[1]
  # which(EventsFound$EventType=="Alternative First Exon")[1]
  # which(EventsFound$EventType=="Alternative Last Exon")[1]
  # which(EventsFound$EventType=="Complex Alt 5' Splice Site")[1]
  # which(EventsFound$EventType=="Complex Alt 3' Splice Site")[1]
  
  try(
    GRseq3 <- foreach(i = seq_len(nrow(EventsFound)),.packages = "GenomicRanges") %dopar%
      # suppressWarnings(GRseq <- foreach(i = seq_len(12),.combine = "c") %dopar%  
      {
        
        # i <- 1  #retained intron
        # i <- 31 #cassette exon
        # i <- 23 #A3'
        # i <- 7 #A5'
        # i <- 4402 #A5' con strand negativo
        # i <- 2 #complex 5' and cassette exon
        # i <- 380 #mutually exclusive
        # i <- 4 #alternative first
        # i <- 27 #alternative last
        # i <- 5 #complex 5'
        # i <- 13 #complex 3'
        
        mygene <- EventsFound$GeneName[i]
        mygeneid <- EventsFound$GeneID[i]
        mytype <- EventsFound$EventType[i]
        myeventsID <- EventsFound$EventID[i]
        myp1 <- EventsFound$Path.1[i]
        myp2 <- EventsFound$Path.2[i]
        mypref <- EventsFound$Path.Reference[i]
        myChr <- gsub(":.*","",EventsFound$GPos[i])
        jjx <- which(names(SG_List)==mygene)
        myEdges <- SG_List[[jjx]]$Edges
        my_strand <- myEdges$Strand[1]
        
        #retained intron, cassette exon, alternative 3', alternative 5'
        if(mytype %in% typeA){
          X <- strsplit(unlist(strsplit(myp1,",")),"-")
          # X
          u1 <- which(myEdges$From == paste0(X[[1]][1],".b") & myEdges$To == paste0(X[[1]][2],".a"))
          u2 <- which(myEdges$From == paste0(X[[3]][1],".b") & myEdges$To == paste0(X[[3]][2],".a"))
          centers <- c(myEdges$Start[u1],myEdges$End[u1],myEdges$Start[u2],myEdges$End[u2])
          
          # alt 5' #####
          ##A5'+ y A3'-  (tenemos que cambiar las anotaciones??):
          if((mytype == "Alternative 5' Splice Site" & my_strand == '+') | 
             (mytype == "Alternative 3' Splice Site" & my_strand == '-')){
            centers <- centers[-2]
            
            nt_intron_1 <- nt_intron
            nt_exon_2 <- nt_exon
            distances <- diff(centers)
            if(distances[2] < nt_intron_1*2){
              nt_intron_1 <- floor(distances[2]/2)
            }
            if(distances[1] < nt_exon_2*2){
              nt_exon_2 <- floor(distances[1]/2)
            }
            
            
            centers <- rep(centers,each=2)
            newstart <- centers-c(nt_exon,-1,
                                  nt_exon_2,-1,
                                  nt_intron_1,-1)
            newend   <- centers+c(0,nt_exon_2-1,
                                  0,nt_intron_1-1,
                                  0,nt_exon)
            nc <- length(centers)
            GR <- data.frame(seqname = rep(myChr,nc),
                             start = newstart,
                             end = newend,
                             strand = myEdges$Strand[1])
            if(my_strand == '+'){
              GR$type <- c("C","D","E","F","G","H")  
            }else{
              GR$type <- c("A","B","C","D","E","F")[6:1]
            }
            # GR <- GRanges(GR)
            # GR 
          }
          ##
          
          # alt 3' #####
          ##A3'+ y A5'-  (tenemos que cambiar las anotaciones??):
          if((mytype == "Alternative 3' Splice Site" & my_strand == '+') | 
             (mytype == "Alternative 5' Splice Site" & my_strand == '-')){
            # centers <- centers_2
            centers <- centers[-4]
            nt_intron_1 <- nt_intron
            nt_exon_2 <- nt_exon
            distances <- diff(centers)
            if(distances[1] < nt_intron_1*2){
              nt_intron_1 <- floor(distances[1]/2)
            }
            if(distances[2] < nt_exon_2*2){
              nt_exon_2 <- floor(distances[2]/2)
            }
            
            centers <- rep(centers,each=2)
            newstart <- centers-c(nt_exon,-1,
                                  nt_intron_1,-1,
                                  nt_exon_2,-1)
            newend <- centers+c(0,nt_intron_1,
                                0,nt_exon_2-1,
                                0,nt_exon)
            
            nc <- length(centers)
            GR <- data.frame(seqname = rep(myChr,nc),
                             start = newstart,
                             end = newend,
                             strand = myEdges$Strand[1])
            if(my_strand == '+'){
              GR$type <- c("A","B","C","D","E","F")  
            }else{
              GR$type <- c("C","D","E","F","G","H")[6:1]
            }
            
            # GR <- GRanges(GR)
            # GR
          }
          
          ##
          # retained intron #####
          ##retained intron:
          if(mytype == "Retained Intron"){
            distances <- diff(centers)[2]
            nt_exon_2 <- nt_exon
            my_distances <- nt_exon_2
            if(distances < my_distances*2){
              nt_exon_2 <- floor(distances/2)
            }
            
            newstart <- centers-c(nt_exon,0,
                                  nt_exon_2,0)
            newend   <- centers+c(0,nt_exon_2-1,
                                  0,nt_exon)
            nc <- length(centers)
            GR <- data.frame(seqname = rep(myChr,nc),
                             start = newstart,
                             end = newend,
                             strand = myEdges$Strand[1])
            if(my_strand == '+'){
              GR$type <- c("C","D","E","F")  
            }else{
              GR$type <- c("C","D","E","F")[4:1]
            }
            
            # GR <- GRanges(GR)
            # GR
          }
          ##
          
          
          # cassette exon: ####
          if(mytype == "Cassette Exon"){
            # centers <- centers_2
            nt_intron_1 <- nt_intron
            nt_intron_2 <- nt_intron
            nt_exon_2 <- nt_exon
            distances <- diff(centers)
            my_distances <- c(nt_intron_1,nt_exon_2,nt_intron_2)
            if(any(distances < my_distances*2)){
              ddx <- which(distances < my_distances*2)
              for(ggx in 1:length(ddx)){
                my_distances[ddx[ggx]] <-  floor(distances[ddx[ggx]]/2)
              }
              nt_intron_1 <- my_distances[1]
              nt_exon_2 <- my_distances[2]
              nt_intron_2 <- my_distances[3]
            }
            
            
            centers <- rep(centers,each=2)
            centers
            newstart <- centers-c(nt_exon,-1, #ab
                                  nt_intron_1,-1, #cd
                                  nt_exon_2,-1, #ef
                                  nt_intron_2,-1) #gh
            
            newend   <- centers+c(0,nt_intron_1-1, #ab
                                  0,nt_exon_2-1, #cd
                                  0,nt_intron_2-1, #ef
                                  0,nt_exon) #gh
            nc <- length(centers)
            GR <- data.frame(seqname = rep(myChr,nc),
                             start = newstart,
                             end = newend,
                             strand = myEdges$Strand[1])
            if(my_strand == '+'){
              GR$type <- c("A","B","C","D","E","F","G","H")
            }else{
              GR$type <- c("A","B","C","D","E","F","G","H")[8:1]
            }
            # GR <- GRanges(GR)
            # GR
          }
          
          
          
        }
        
        # mutually exclusive exons ####
        if(mytype == "Mutually Exclusive Exons"){
          
          X <- strsplit(unlist(strsplit(myp1,",")),"-")
          u1 <- which(myEdges$From == paste0(X[[1]][1],".b") & myEdges$To == paste0(X[[1]][2],".a"))
          u2 <- which(myEdges$From == paste0(X[[3]][1],".b") & myEdges$To == paste0(X[[3]][2],".a"))
          Y <- strsplit(unlist(strsplit(myp2,",")),"-")
          u3 <- which(myEdges$From == paste0(Y[[2]][1],".a") & myEdges$To == paste0(Y[[2]][2],".b"))
          centers <- c(myEdges$Start[u1],myEdges$End[u1],
                       myEdges$Start[u2],myEdges$End[u2],
                       myEdges$Start[u3],myEdges$End[u3])
          centers <- sort(centers)
          distances <- diff(centers)
          nt_exon_1 <- nt_exon
          nt_exon_2 <- nt_exon
          nt_intron_1 <- nt_intron
          nt_intron_2 <- nt_intron
          nt_intron_3 <- nt_intron
          my_distances <- c(nt_intron_1,nt_exon_1,nt_intron_2,nt_exon_2,nt_intron_3)
          if(any(distances < my_distances*2)){
            ddx <- which(distances < my_distances*2)
            for(ggx in 1:length(ddx)){
              my_distances[ddx[ggx]] <-  floor(distances[ddx[ggx]]/2)
            }
            nt_intron_1 <- my_distances[1]
            nt_exon_1 <- my_distances[2]
            nt_intron_2 <- my_distances[3]
            nt_exon_2 <- my_distances[4]
            nt_intron_3 <- my_distances[5]
          }
          
          centers <- rep(centers,each=2)
          centers
          newstart <- centers-c(nt_exon,-1,
                                nt_intron_1,-1,
                                nt_exon_1,-1,
                                nt_intron_2,-1,
                                nt_exon_2,-1,
                                nt_intron_3,-1)
          newend   <- centers+c(0,nt_intron_1-1,
                                0,nt_exon_1-1,
                                0,nt_intron_2-1,
                                0,nt_exon_2-1,
                                0,nt_intron_3-1,
                                0,nt_exon-1)
          nc <- length(centers)
          GR <- data.frame(seqname = rep(myChr,nc),
                           start = newstart,
                           end = newend,
                           strand = myEdges$Strand[1])
          if(my_strand == '+'){
            GR$type <- c("A","B","C","D","E","F","C","D","E","F","G","H")
          }else{
            GR$type <- c("A","B","C","D","E","F","C","D","E","F","G","H")[12:1]
          }
          # GR <- GRanges(GR)
          # GR
        }
        
        # alternative first exon ####
        if((mytype == "Alternative First Exon" & my_strand == "+") |
           (mytype == "Alternative Last Exon" & my_strand == "-")){
          myp1
          myp2
          mypref
          X <- strsplit(unlist(strsplit(myp1,",")),"-")
          Y <- strsplit(unlist(strsplit(myp2,",")),"-")
          X <- do.call(rbind,X)
          Y <- do.call(rbind,Y)
          
          i1 <- which(X[,2] %in% Y[,2])
          i2 <- which(Y[,2] %in% X[,2])
          
          if(length(i1)==0 | length(i2)==0){
            return(NULL)
          }
          
          u1 <- c(X[i1,1],Y[i2,1],X[i1,2])
          p1 <- which(myEdges$From==paste0(u1[1],".a"))
          p2 <- which(myEdges$From==paste0(u1[2],".a"))
          pref <- which(myEdges$From==paste0(u1[3],".a"))
          centers <- c(myEdges[p1,c(4,5)],myEdges[p2,c(4,5)],myEdges[pref,c(4)])
          centers <- as.vector(unlist(centers))
          centers <- sort(centers)
          
          distances <- diff(centers)
          nt_exon_1 <- nt_exon
          nt_exon_2 <- nt_exon
          nt_intron_1 <- nt_intron
          nt_intron_2 <- nt_intron
          my_distances <- c(nt_exon_1,nt_intron_1,nt_exon_2,nt_intron_2)
          if(any(distances < my_distances*2)){
            ddx <- which(distances < my_distances*2)
            for(ggx in 1:length(ddx)){
              my_distances[ddx[ggx]] <-  floor(distances[ddx[ggx]]/2)
            }
            nt_exon_1 <- my_distances[1]
            nt_intron_1 <- my_distances[2]
            nt_exon_2 <- my_distances[3]
            nt_intron_2 <- my_distances[4]
          }
          
          
          
          centers <- rep(centers,each=2)
          centers
          newstart <- centers-c(nt_intron,-1,
                                nt_exon_1,-1,
                                nt_intron_1,-1,
                                nt_exon_2,-1,
                                nt_intron_2,-1)
          newend   <- centers+c(0,nt_exon_1-1,
                                0,nt_intron_1-1,
                                0,nt_exon_2-1,
                                0,nt_intron_2-1,
                                0,nt_exon)
          nc <- length(centers)
          GR <- data.frame(seqname = rep(myChr,nc),
                           start = newstart,
                           end = newend,
                           strand = myEdges$Strand[1])
          if(my_strand == '+'){
            GR$type <- c("C","D","E","F","C","D","E","F","G","H")
          }else{
            GR$type <- c("A","B","C","D","E","F","C","D","E","F")[10:1]
          }
          # GR <- GRanges(GR)
          # GR
        }
        
        # alternative last exons ####
        if((mytype == "Alternative Last Exon" & my_strand == "+") |
           (mytype == "Alternative First Exon" & my_strand == "-")){
          myp1
          myp2
          mypref
          X <- strsplit(unlist(strsplit(myp1,",")),"-")
          Y <- strsplit(unlist(strsplit(myp2,",")),"-")
          X <- do.call(rbind,X)
          Y <- do.call(rbind,Y)
          
          i1 <- which(X[,1] %in% Y[,1])
          i2 <- which(Y[,1] %in% X[,1])
          
          if(length(i1)==0 | length(i2)==0){
            return(NULL)
          }
          
          u1 <- c(X[i1,1],Y[i2,2],X[i1,2])
          p1 <- which(myEdges$From==paste0(u1[3],".a"))
          p2 <- which(myEdges$From==paste0(u1[2],".a"))
          pref <- which(myEdges$From==paste0(u1[1],".a"))
          centers <- c(myEdges[pref,c(5)],myEdges[p1,c(4,5)],myEdges[p2,c(4,5)])
          centers <- as.vector(unlist(centers))
          
          centers <- sort(centers)
          distances <- diff(centers)
          nt_exon_1 <- nt_exon
          nt_exon_2 <- nt_exon
          nt_intron_1 <- nt_intron
          nt_intron_2 <- nt_intron
          my_distances <- c(nt_intron_1,nt_exon_1,nt_intron_2,nt_exon_2)
          if(any(distances < my_distances*2)){
            ddx <- which(distances < my_distances*2)
            for(ggx in 1:length(ddx)){
              my_distances[ddx[ggx]] <-  floor(distances[ddx[ggx]]/2)
            }
            nt_intron_1 <- my_distances[1]
            nt_exon_1 <- my_distances[2]
            nt_intron_2 <- my_distances[3]
            nt_exon_2 <- my_distances[4]
          }
          
          centers <- rep(centers,each=2)
          # centers
          newstart <- centers-c(nt_exon,-1,
                                nt_intron_1,-1,
                                nt_exon_1,-1,
                                nt_intron_2,-1,
                                nt_exon_2,-1)
          newend   <- centers+c(0,nt_intron_1-1,
                                0,nt_exon_1-1,
                                0,nt_intron_2-1,
                                0,nt_exon_2-1,
                                0,nt_intron-1)
          nc <- length(centers)
          GR <- data.frame(seqname = rep(myChr,nc),
                           start = newstart,
                           end = newend,
                           strand = myEdges$Strand[1])
          if(my_strand == '+'){
            GR$type <- c("A","B","C","D","E","F","C","D","E","F")
          }else{
            GR$type <- c("C","D","E","F","C","D","E","F","G","H")[10:1]
          }
          # GR <- GRanges(GR)
          # GR
        }
        
        
        # GR2 <- reduce(GR)
        h <- nrow(GR)
        GR$SeqID <- paste0(myeventsID,"-",1:h)
        GR$EventID <- myeventsID
        GR$Length.Seq <- GR$end-GR$start
        if(any(GR$Length.Seq < 3)){
          rrx <- which(GR$Length.Seq < 1)
          GR <- GR[-rrx,,drop=F]
        }
        
        return(GR)
      }
  )
  stopCluster(cl)
  GRseq3 <- Filter(Negate(is.null), GRseq3)
  GRseq3 <- do.call(rbind,GRseq3)
  head(GRseq3)
  GRseq3 <- GRanges(seqnames = GRseq3$seqname,
                    ranges = IRanges(GRseq3$start,GRseq3$end),
                    strand = GRseq3$strand,
                    SeqID=GRseq3$SeqID,
                    EventID=GRseq3$EventID,
                    Length.Seq=GRseq3$Length.Seq,
                    type=GRseq3$type)
  
  return(GRseq3)
  
}









