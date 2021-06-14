#!/usr/bin/env Rscript
##############################################################################
##                           Description                                   ##
##############################################################################
#Title: Heat map of V an J segments
#Description:comparation B6 and NOD 
#Installing and loading required packages
if(!require(gplots))
{install.packages("gplots")
}
if(!require(RColorBrewer))
{install.packages("RColorBrewer")
}
if(!require(ggplot2))
{install.packages("ggplot2")
}
if(!require(plotrix))
{install.packages("plotrix")
}
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(plotrix)
#################################
#Take arguments from command line:
args = commandArgs(trailingOnly=TRUE)
##############################################################################
##                           Initial Setup                                 ##
##############################################################################
#set working directories
#setwd(args[1])
Sample_names <- c("B6", "NOD")
#READ file:
df_V<-read.table(args[2], header=T, sep="\t", comment="")
df_J<-read.table(args[3], header=T, sep="\t", comment="")
#Change names:
df_V$sample_id <- Sample_names
df_J$sample_id <- Sample_names
##############################################################################
##                          Process data                                    ##
##############################################################################
#Vsegments
v_col_count <- as.numeric(ncol(df_V)-2)
vcols<-(ncol(df_V)-v_col_count + 1):ncol(df_V)
df_V[, vcols] <- apply(df_V[, vcols], 2, as.numeric)
V <- as.matrix(df_V[, vcols])
#J segments
j_col_count <- as.numeric(ncol(df_J)-2)
jcols<-(ncol(df_J)-j_col_count + 1):ncol(df_J)
df_J[, jcols] <- apply(df_J[, jcols], 2, as.numeric)
J <- as.matrix(df_J[, jcols])
##############################################################################
##                          Graphics                                        ##
##############################################################################
#Vsegments
svg(args[4])
heatmap.2(t(V),  na.rm=F, labCol = df_V[,"sample_id"],
          na.col="grey50",
          col=colorRampPalette(c("#2f98ce", "#e0f3db", "#f47104")),
          density.info="none", trace="none",scale="column", srtCol=0, cexCol = 1.2)
dev.off()

#Vsegments
svg(args[5])
heatmap.2(t(J),  na.rm=F, labCol = df_J[,"sample_id"],
          na.col="grey50",
          col=colorRampPalette(c("#2f98ce", "#e0f3db", "#f47104")),
          density.info="none", trace="none",scale="column", srtCol=0, cexCol = 1.2)
dev.off()