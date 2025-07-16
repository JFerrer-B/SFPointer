

source("./SF_pointer_functions/Create_ExS_matrix.r")

#firs step: create EventRegions (done) #######

#to run the followig code you need the files returned by eventpointer detection function.

# tablainfoevents <- "path files /EventsFound_gencode24.txt"
# EventsFound <- read.delim(file=tablainfoevents,stringsAsFactors = FALSE)
# EventsFound$EventID <- paste0(EventsFound$GeneName,"_",EventsFound$EventNumber)
# 
# load("path files/EventsXtrans_gc24.RData")
# SG_List <- EventsXtrans_gc24$SG_List
# 
# Events_Regions <- Compute_Events_Regions(EventsFound=EventsFound_2,
#                                          SG_List=SG_List,
#                                          cores=1,
#                                          nt_intron=300,nt_exon=100)


#note that this functions depends on: the files returned by EventPointer and nt_inron and nt_exon.
#if you don't change this parameters this function you need to run this function only once.


#second step: create ExS matrix ############

load("./data/data_create_ExS/Events_Regions.RData")
load("./data/data_create_ExS/peaks_fus_mbnl1_example.RData")

ExS_test <- CreateExSmatrix(POSTAR = peaks,
                            EventsRegions = Events_Regions)
head(ExS_test)
table(ExS_test[,1])
table(ExS_test[,2])







