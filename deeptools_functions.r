
require("dplyr")
require("ggplot2")
require("reshape2")
require("rtracklayer")
library("BSgenome.Hsapiens.UCSC.hg19")
library("plyranges")
library("tidyverse")
library("RColorBrewer")
seqlens = seqlengths( Hsapiens );

##################
## Functions for plotting deeptools
#################

## PATH OF WORKING DIRECTORY
working_dir <- "/home/scollins/Heatmaps_and_profiles/"

## PATH TO `COMPUTE MATRIX` FUNCTION OF DEEPTOOLS
deeptoolspath <- "/home/scollins/anaconda3/envs/deeptools/bin/computeMatrix"

## VOLUMES TO MOUNT WITH RSHINY APP 
volumes <- c(DATA_HOME='/home/scollins/')


## FUNCTION TO READ DEEPTOOLS MATRIX AND TRANSFORM TO GGPLOT DATAFRAME
plot_group <- function(FILES, vars, regions, strand="same"){
  
  dat.plot <- NULL
  dat.heatmap <- NULL
  dat.boxplot <- NULL
  
  
  message("Reading in data...")
  ## load in data to dfs
  for(n in 1:length(FILES)){
    df <- read.csv(FILES[n], sep="\t", header = FALSE, skip=1)
    rownames(df) <- df$V4
    #rownames(df) <- paste0(df$V1,":",df$V2,"-",df$V3)
    df <- (df[,-c(1:6)])
    df[df=="NaN"] <- 0
    
    if(strand=="split"){
      df <- df[,ncol(df):1]
    }else{
      df <- df
    }
    
      if(nrow(df)>1000){
        s <- nrow(df) / 1000
        df <- df[ order(rowMeans(df)), ]
        df <- aggregate(df, list(rep(1:(nrow(df) %/% s + 1), each = s, len = nrow(df))), mean)[-1];
      }
    
    sub.dat.plot <- data.frame(
      Window=c(1:ncol(df)),
      Value=colMeans((df)),
      variable = vars[n],
      Peaks=regions[n]
    )
    
    if(is.null(dat.plot)){
      dat.plot <- sub.dat.plot
    }else{
      dat.plot <- rbind(dat.plot,sub.dat.plot)
    }
    
    
    colnames(df) <- 1:ncol(df)
    melted <- melt(df)
    colnames(melted) <- c("window","value")
    melted$variable <- vars[n]
    melted$Peaks <- regions[n]
    melted$site <- rownames(df)
    
    if(is.null(dat.heatmap)){
      dat.heatmap <- melted
    }else{
      dat.heatmap <- rbind(dat.heatmap,melted)
    }   
  
    
  sub.dat.boxplot <- data.frame(
    Value=rowMeans(df),
    variable = vars[n],
    Peaks=regions[n],
    site=rownames(df)
  )
  
  if(is.null(dat.plot)){
    dat.boxplot <- sub.dat.boxplot
  }else{
    dat.boxplot <- rbind(dat.boxplot,sub.dat.boxplot)
  }
  
  }
  
  return(list(dat.heatmap,dat.plot,dat.boxplot))
  
} 

