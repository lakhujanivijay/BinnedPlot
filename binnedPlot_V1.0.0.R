###############
##Author : Chirag Parsania
##email : chirag.parsania@gmail.com /yb57653@umac.mo
###############

#############################################################
## Require 2 files to run this script 
## ----------------------------------------------------
## file 1 genome file : file containing genomic varibale defined along genome E.g coverage or depth 
## file 2 feature file : any genomic features E.g genes/mrna/cds etc. 
############################################################
## file 1 format
##----------------------------------------------------
## must have following columns without headers. 
## chr start end score 
## start is always smaller then end and there must not be any overlap between start and end.
## check overlap using isDisjoint command.
############################################################
## file 2 format
##----------------------------------------------------
## It's a feature file, deleminated by tab.
## must have folllwing columns without headers. 
## chr	start	end	strand	score	name
############################################################



#############################################
## load function 
#############################################

getBinnedMatrix <- function(genomeVaribaleFile, featureFile , flankingBp = 500 , numberOfBinsInGeneBody = 100 , numberOfBinsInFlanking = 25 ,isVariableOVerlapped = TRUE, isFeatureDF = FALSE,isGenomeRLEobject = TRUE){
        ## required libraries 
        
        library(GenomicRanges)
        library(IRanges)
        
        ##########################################
        ## genomeVaribaleFile = .bdg file or RLE object.  isGenomeRLEobject = TRUE if RLE is provided. 
        ## featureFile = feature file or dataframe of feature file. isFeatureDF = TRUE if feature dataframe provided. 
        ## flankingBp = integer number, default 500
        ## numberOfBinsInGeneBody = integer number, default 100
        ## numberOfBinsInFlanking = integer number, default 25
        ## isVariableOVerlapped = boolean, default TRUE. TRUE indicates genome file has one base pair overlap. So increase start and end by one bp to remove overlap.
        ## isFeatureDF = boolean, default FALSE. FALSE indicates your features are in file object and not in dataframe.
        ## isGenomeRLEobject = boolean, default FALSE. FALSE indicates genome provided as .bdg file 
        
        
        #############################################
        ## prepare genome data --> 
        #############################################
        
        
        #############################################
        ## isGenomeRLEobject : FALSE -->  create RLE object from given bdg file. 
        ## if user have same bdg file and  multiple feature lists then evry time it requires same RLE object from .bdg file 
        ## Create  RLE object once and next time directly give RLE object rather creating it again. 
        ## isGenomeRLEobject : TRUE --> directly use RLE provided RLE object  
        #############################################
        
        ###################Test Varibale########################
        # genomeVaribaleFile <- "/Users/chiragparsania/Documents/Projects/1_My_PhD/Joshua_Lab/Sc_polII_data_pooja/Sc_YPD_Pol2.counted.bdg.mapped"
        # genomeVaribaleFile <-"/Users/chiragparsania/Documents/Projects/Shuhui/A_nidulans_FGSC_A4_version_s10-m03-r13.fasta_bendabilty.bdg.1"
        # genomeVaribaleFile <- genomeVariableFiles
        # 
        #featureFile <- "testFeatures.txt"
        # flankingBp = 500
        # numberOfBinsInGeneBody = 200
        # numberOfBinsInFlanking = 50
        # isVariableOVerlapped = FALSE
        # isFeatureDF = FALSE
        # isGenomeRLEobject = FALSE
        # # ######################################################
        
        ############################################
        if(!isGenomeRLEobject){
                print("Create RLE object from BDG file. It will take a moment...")
                genome_bdg <- read.table(genomeVaribaleFile,sep="\t",quote = "")
                genome_bdg <- genome_bdg[,c(1:4)]
                colnames(genome_bdg) <- c("chr","start","end","score")
                
                ## increase  start/end position by 1 incase overlapping regions. 
                
                if(isVariableOVerlapped){
                        genome_bdg$start <- genome_bdg$start +1
                        #genome_bdg$end <- genome_bdg$end + 1
                }
                
                varibale_ranges <- makeGRangesFromDataFrame(genome_bdg,keep.extra.columns = T)
                
                varibale_rle  <- mcolAsRleList(varibale_ranges ,"score")
                
        }else{
                print("RLE object provided. It will make process faster !")
                varibale_rle <- genomeVaribaleFile
        }
        #############################################
        ## prepare feature data : data from file or dataframe object.
        #############################################
        if(!isFeatureDF){
                genes <- read.table(featureFile,sep="\t",quote = "",header = T)
                head(genes)
                genes_ranges <- makeGRangesFromDataFrame(genes,keep.extra.columns = T)
        }else{
                genes_ranges <- makeGRangesFromDataFrame(featureFile,keep.extra.columns = T)   
        }
        #############################################
        ## Do tile of each feature 
        #############################################
        
        binnedFeatures <- lapply(genes_ranges , function(i){
                totalLength <- (flankingBp+width(i)+flankingBp)
                totalBins  <- (numberOfBinsInFlanking + numberOfBinsInGeneBody + numberOfBinsInFlanking)
                
                ### start and end for flanking regions : start is always smaller then end. 
                
                geneStart = min(start(i),end(i))
                geneEnd = max(start(i),end(i))
                
                ## gene body bins 
                
                geneBodyRanges = IRanges(start = seq(geneStart,to = geneEnd , length.out = numberOfBinsInGeneBody),
                                         width = width(i)/numberOfBinsInGeneBody)
                
                ## flank 1 
                flank1 = IRanges(start = seq(from=(geneStart - flankingBp) ,to = geneStart , length.out = numberOfBinsInFlanking),
                                 width = flankingBp/numberOfBinsInFlanking)
                ## flank 2 
                flank2 = IRanges(start = seq(from=geneEnd ,to = (geneEnd + flankingBp) , length.out = numberOfBinsInFlanking),
                                 width = flankingBp/numberOfBinsInFlanking)
                
                ## get starts of all ranges.
                starts <- c(start(flank1),start(geneBodyRanges),start(flank2))
                
                gr <- GRanges(seqnames = Rle(rep(as.character(seqnames(i)),totalBins)),
                              strand = Rle(rep(as.character(strand(i)),totalBins)),
                              ranges = IRanges(start = starts,
                                               width = totalLength/totalBins))
                
                seqlevels(gr) <- names(varibale_rle)
                tt <- binnedAverage(gr, varibale_rle, "binned_var")
                mcols(tt)$name <- rep(mcols(i)$name, length(tt))
                return(tt)
        })
        
        head(binnedFeatures)
        
        ## calculate binned matrix. if strand is - , reverse the orientation.
        
        binnedMat <- sapply(binnedFeatures, function(features){
                ## get strand of first range. Each elemnt of list contain ranges for given feature. So all ranges must have same strand.
                range1 <-  features[1]
                strnd <- as.character(strand(range1))
                binnedval <- mcols(features)$binned_var
                gName  <- as.character(mcols(features)$name[1])
                ## if strand is negative then print it reverse 
                if(strnd == "-"){
                        return(rev(binnedval))
                }else{
                        return(binnedval)
                }
        })
        binnedMat <- t(binnedMat)
        return(binnedMat)
}
