########################################################
## This script explain how to run "binnedPlot.R" function. 
########################################################

########################################################
## "binnedPlot.R" can be run for 
##      1). One featureset One bedGraph file (one experiment) 
##      2). One featureset multiple bedGraph file (multiple experiment)
##      3). multiple featureset multiple bedGraph file 
########################################################

########################################################
## 1). one featureset one bedgraph file 
########################################################

### set path of binnedPlot.R here. Default is current directory 

source("./binnedPlot.R")
source("./plotMultipleGplots.R")
library(ggplot2)
library(gplots)
library(GenomicRanges)
library(reshape2)
#########################
## Prepare sample data 
#########################

########################################################
## 1) bed graph file 
## Given bed graph files are taken from published data 
## segal et. al Nature 458, 362-366 (19 March 2009) | doi:10.1038/nature07667
########################################################

genomeVaribaleFilesPath <- "./sampleData/bdg"
genomeVariableFiles<-list.files(genomeVaribaleFilesPath,pattern = "*.bdg$",full.names = T)

########################################################
## 1) feature file  
## given feature file contains total number of s.cerevisiae genes. 
## approx. 6k genes devided in to two cluster by column name "clust"
########################################################

## prepare gene features. Incase of multiple featureset group them by colname "clust"

featureFile = "./sampleData/featureSet/sampleFeatures.txt"
featureDf <- read.table(featureFile,sep="\t",quote = "",header = T)

## prepare list  if multiple featureSet.
featureAsList <- split(featureDf,f = featureDf$clust) 
print(c("number of feature set are",length(featureAsList)),quote =F)

#featureAsList <- list(featureDf)

### create nuclesome RLE objects. Creating RLE obkect first is useful if you have multiple bedgraph files.

# Start the clock!
ptm <- proc.time()

genome_varibale_rle_objects <- lapply(genomeVariableFiles, function(bdg){
        print(paste("Reading",basename(bdg),"...",sep = " "),quote = F,print.gap = 0)
        
        genome_bdg <- read.table(bdg,sep="\t",quote = "")
        colnames(genome_bdg) <- c("chr","start","end","score")
        varibale_ranges <- makeGRangesFromDataFrame(genome_bdg,keep.extra.columns = T)
        
        if(!(isDisjoint(varibale_ranges))){
                genome_bdg$start <- genome_bdg$start+1
                varibale_ranges <- makeGRangesFromDataFrame(genome_bdg,keep.extra.columns = T)
        }
        varibale_rle  <- mcolAsRleList(varibale_ranges ,"score")
        return(varibale_rle)
})

# Stop the clock
ptm <- proc.time()

gplots = NULL
for (i in c(1:length(featureAsList))){
        # start the clock 
        ptm <- proc.time()
        
        selectedFeatureClust <- featureAsList[[i]]
        
        ## for each NOD (nuclesome occupancy data) generate binned matrix 
        
        print(paste("calculating binned matrix for ", names(featureAsList)[i]))
        
        ## supply list of genomeVaribale File. For genomeVaribale binnedMatrix List will be created
        binnedMat <- lapply(genome_varibale_rle_objects, function(X){
                return(getBinnedMatrix(genomeVaribaleFile = X, featureFile = selectedFeatureClust,numberOfBinsInGeneBody = 200,numberOfBinsInFlanking = 50,flankingBp = 500,isVariableOVerlapped = FALSE,isFeatureDF = TRUE,isGenomeRLEobject=TRUE))        
        })
        
        ## remove genes which have NA values 
        
        naRemoved <- lapply(binnedMat, function(dat){
                isNA <- apply(dat,1,function(X){any(is.na(X))})
                geneOfInterest <- dat[!isNA,]
                return(geneOfInterest)
        })
        
        ## get col mean of each dataset.
        
        avgAtEachPos <- lapply(naRemoved, function(dat){
                return(colMeans(dat))
        }) 
        ## prepare data fro ggplot 
        
        names(avgAtEachPos) <- basename(genomeVariableFiles)
        melted<- melt(avgAtEachPos)
        
        ## for each avaraged variable create index attached to each genome varible file. 
        chk <- cbind(melted,position=rep(1:300, length(genomeVariableFiles)))
        
        print("preparing plot")
        
        p2 = ggplot(chk, aes(x=position,y=value, group = L1, colour=L1)) + geom_line()+
                theme(axis.text.x= element_text(face="bold", colour="brown", size=14),
                      axis.text.y = element_text(face="bold", color="#993333",size=14, angle=45))+
                labs(y="Counts")+
                ggtitle(paste("clust" , i , sep=""))+theme_bw()
        gplots[[i]] <- p2
        elapsed <- (proc.time() - ptm)["elapsed"]
        print(paste("processing time for ",featureAsList[[i]],"is",elapsed,sep = " "))
        
}    

## plot data 

multiplot(plotlist = gplots, cols = 2)

## process raw gplots. Change plot parametrs as per your requirements

minYlim <- -2
maxYlim <- 0.5

gplots_with_ylim <- lapply(gplots , function(plt){return(plt + coord_cartesian(ylim=c(minYlim, maxYlim)))})
gplots_with_label_and_title <- lapply(gplots_with_ylim , function(plt){return( plt + scale_x_continuous(breaks = seq(from=0,to=300,by=50),labels=c("-100","TSS","","","","TES","+100")) 
                                                                               + ylab(expression(atop("Average ",paste("Nulcesome occupancy")))))})
multiplot(plotlist = gplots_with_label_and_title, cols = 2)

