#!/usr/bin/env Rscript
##############################################################################
##                           Description                                   ##
##############################################################################
#Title: Reads, UMIs and Clonotypes Histograms
#Description: 
#Installing and loading required packages
if(!require(stringr))
{install.packages("stringr")
}
library(stringr)
#################################
#Take arguments from command line:
args = commandArgs(trailingOnly=TRUE)
##############################################################################
##                           Initial Setup                                 ##
##############################################################################
#set working directories
#setwd(working_directory)
#READ file:
UvR1_Cond1 <- read.table(args[3], header = FALSE)
UvR2_Cond1 <- read.table(args[4], header = FALSE)
UvR1_Cond2 <- read.table(args[5], header = FALSE)
UvR2_Cond2 <- read.table(args[6], header = FALSE)
UvC1_Cond1 <- read.table(args[7], header = TRUE)
UvC2_Cond1 <- read.table(args[8], header = TRUE)
UvC1_Cond2 <- read.table(args[9], header = TRUE)
UvC2_Cond2 <- read.table(args[10], header = TRUE)
###############################################################################
##                          Process Data ALL SAMPLES                         ##
###############################################################################
UvC1_Cond1 <- UvC1_Cond1[1]
UvC2_Cond1 <- UvC2_Cond1[1]
UvC1_Cond2 <- UvC1_Cond2[1]
UvC2_Cond2 <- UvC2_Cond2[1]
colnames(UvR1_Cond1) <- c("count") 
colnames(UvR2_Cond1) <- c("count")
colnames(UvR1_Cond2) <- c("count") 
colnames(UvR2_Cond2) <- c("count")
#Add sequential numbers:
UvC1_Cond1$seq <- seq(nrow(UvC1_Cond1), 1)
UvC2_Cond1$seq <- seq(nrow(UvC2_Cond1), 1)
UvC1_Cond2$seq <- seq(nrow(UvC1_Cond2), 1)
UvC2_Cond2$seq <- seq(nrow(UvC2_Cond2), 1)
UvR1_Cond1$seq <- seq(nrow(UvR1_Cond1), 1)
UvR2_Cond1$seq <- seq(nrow(UvR2_Cond1), 1)
UvR1_Cond2$seq <- seq(nrow(UvR1_Cond2), 1) 
UvR2_Cond2$seq <- seq(nrow(UvR2_Cond2), 1)
### Percentages of minimum values
min_UvC1_Cond1 <- round((sum(UvC1_Cond1$count == 1)*100)/nrow(UvC1_Cond1), digits = 0)
min_UvC2_Cond1 <- round((sum(UvC2_Cond1$count == 1)*100)/nrow(UvC2_Cond1), digits = 0)
min_UvC1_Cond2 <- round((sum(UvC1_Cond2$count == 1)*100)/nrow(UvC1_Cond2), digits = 0)
min_UvC2_Cond2 <- round((sum(UvC2_Cond2$count == 1)*100)/nrow(UvC2_Cond2), digits = 0)
min_UvR1_Cond1 <- round((sum(UvR1_Cond1$count == 1)*100)/nrow(UvR1_Cond1), digits = 0)
min_UvR2_Cond1 <- round((sum(UvR2_Cond1$count == 2)*100)/nrow(UvR2_Cond1), digits = 0)
min_UvR1_Cond2 <- round((sum(UvR1_Cond2$count == 1)*100)/nrow(UvR1_Cond2), digits = 0)
min_UvR2_Cond2 <- round((sum(UvR2_Cond2$count == 2)*100)/nrow(UvR2_Cond2), digits = 0)
#Limits for Y axis:
UvC1_Ylim <- as.numeric(max(max(UvC1_Cond1$count), max(UvC1_Cond2$count))) 
UvC2_Ylim <- as.numeric(max(max(UvC2_Cond1$count), max(UvC2_Cond2$count))) 
UvR1_Ylim <- as.numeric(max(max(UvR1_Cond1$count), max(UvR1_Cond2$count))) 
UvR2_Ylim <- as.numeric(max(max(UvR2_Cond1$count), max(UvR2_Cond2$count))) 
###############################################################################
##                                   Graphics                                ##
###############################################################################
### UMI vs Clonotypes:
###UvC1_Cond1
xlimUvC1_Cond1 <- rev(range(UvC1_Cond1$seq))
xlimUvC1_Cond1[1] <- xlimUvC1_Cond1[1]-((xlimUvC1_Cond1[1])*0.03)
xlimUvC1_Cond1[2] <- xlimUvC1_Cond1[2]+((xlimUvC1_Cond1[1])*0.03)
svg(paste(args[1],"_UMIvsClonotype_1read",".svg", sep = ""), width=7, height=7)  
plot(UvC1_Cond1$seq,UvC1_Cond1$count, type="l", xlim = xlimUvC1_Cond1, ylim=c(1.5,UvC1_Ylim), col="dark blue", lwd=7, xlab="Ranked Clonotypes", ylab="Number of MIGS", main=args[1])
text((xlimUvC1_Cond1[1]/5),UvC1_Ylim, paste("Clonotypes with 1 MIG:", min_UvC1_Cond1, "%"))
dev.off()
###UvC2_Cond1
xlimUvC2_Cond1 <- rev(range(UvC2_Cond1$seq))
xlimUvC2_Cond1[1] <- xlimUvC2_Cond1[1]-((xlimUvC2_Cond1[1])*0.03)
xlimUvC2_Cond1[2] <- xlimUvC2_Cond1[2]+((xlimUvC2_Cond1[1])*0.03)
svg(paste(args[1],"_UMIvsClonotype_2read",".svg", sep = ""), width=7, height=7)  
plot(UvC2_Cond1$seq,UvC2_Cond1$count, type="l", xlim = xlimUvC2_Cond1, ylim=c(1.5,UvC2_Ylim), col="dark blue", lwd=7, xlab="Ranked Clonotypes", ylab="Number of MIGS", main=args[1])
text((xlimUvC2_Cond1[1]/5),UvC2_Ylim, paste("Clonotypes with 1 MIG:", min_UvC2_Cond1, "%"))
dev.off()
###UvC1_Cond2
xlimUvC1_Cond2 <- rev(range(UvC1_Cond2$seq))
xlimUvC1_Cond2[1] <- xlimUvC1_Cond2[1]-((xlimUvC1_Cond2[1])*0.03)
xlimUvC1_Cond2[2] <- xlimUvC1_Cond2[2]+((xlimUvC1_Cond2[1])*0.03)
svg(paste(args[2],"_UMIvsClonotype_1read",".svg", sep = ""), width=7, height=7)  
plot(UvC1_Cond2$seq,UvC1_Cond2$count, type="l", xlim = xlimUvC1_Cond2, ylim=c(1.5,UvC1_Ylim), col="dark blue", lwd=7, xlab="Ranked Clonotypes", ylab="Number of MIGS", main=args[2])
text((xlimUvC1_Cond2[1]/5),UvC1_Ylim, paste("Clonotypes with 1 MIG:", min_UvC1_Cond2, "%"))
dev.off()
###UvC2_Cond2
xlimUvC2_Cond2 <- rev(range(UvC2_Cond2$seq))
xlimUvC2_Cond2[1] <- xlimUvC2_Cond2[1]-((xlimUvC2_Cond2[1])*0.03)
xlimUvC2_Cond2[2] <- xlimUvC2_Cond2[2]+((xlimUvC2_Cond2[1])*0.03)
svg(paste(args[2],"_UMIvsClonotype_2read",".svg", sep = ""), width=7, height=7)  
plot(UvC2_Cond2$seq,UvC2_Cond2$count, type="l", xlim = xlimUvC2_Cond2, ylim=c(1.5,UvC2_Ylim), col="dark blue", lwd=7, xlab="Ranked Clonotypes", ylab="Number of MIGS", main=args[2])
text((xlimUvC2_Cond2[1]/5),UvC2_Ylim, paste("Clonotypes with 1 MIG:", min_UvC2_Cond2, "%"))
dev.off()
### UMI vs Reads:
###UvR1_Cond1
xlimUvR1_Cond1 <- rev(range(UvR1_Cond1$seq))
xlimUvR1_Cond1[1] <- xlimUvR1_Cond1[1]-((xlimUvR1_Cond1[1])*0.03)
xlimUvR1_Cond1[2] <- xlimUvR1_Cond1[2]+((xlimUvR1_Cond1[1])*0.03)
svg(paste(args[1],"_UMIvsReads_1read",".svg", sep = ""), width=7, height=7)  
plot(UvR1_Cond1$seq,UvR1_Cond1$count, log = "y", type="l", xlim = xlimUvR1_Cond1, ylim=c(1,UvR1_Ylim), col="dark blue", lwd=7, xlab="Ranked MIG", ylab="Log 10 (Number of Reads)", main=args[1])
text((xlimUvR1_Cond1[1]/5),UvR1_Ylim, paste("MIGs with 1 read:", min_UvR1_Cond1, "%"))
dev.off()
###UvR2_Cond1
xlimUvR2_Cond1 <- rev(range(UvR2_Cond1$seq))
xlimUvR2_Cond1[1] <- xlimUvR2_Cond1[1]-((xlimUvR2_Cond1[1])*0.03)
xlimUvR2_Cond1[2] <- xlimUvR2_Cond1[2]+((xlimUvR2_Cond1[1])*0.03)
svg(paste(args[1],"_UMIvsReads_2read",".svg", sep = ""), width=7, height=7)  
plot(UvR2_Cond1$seq,UvR2_Cond1$count, log = "y", type="l", xlim = xlimUvR2_Cond1, ylim=c(1.5,UvR2_Ylim), col="dark blue", lwd=7, xlab="Ranked MIG", ylab="Log 10 (Number of Reads)", main=args[1])
text((xlimUvR2_Cond1[1]/5),UvR2_Ylim, paste("MIGs with 2 reads:", min_UvR2_Cond1, "%"))
dev.off()
###UvR1_Cond2
xlimUvR1_Cond2 <- rev(range(UvR1_Cond2$seq))
xlimUvR1_Cond2[1] <- xlimUvR1_Cond2[1]-((xlimUvR1_Cond2[1])*0.03)
xlimUvR1_Cond2[2] <- xlimUvR1_Cond2[2]+((xlimUvR1_Cond2[1])*0.03)
svg(paste(args[2],"_UMIvsReads_1read",".svg", sep = ""), width=7, height=7)  
plot(UvR1_Cond2$seq,UvR1_Cond2$count, log = "y", type="l", xlim = xlimUvR1_Cond2, ylim=c(1,UvR1_Ylim), col="dark blue", lwd=7, xlab="Ranked MIG", ylab="Log 10 (Number of Reads)", main=args[2])
text((xlimUvR1_Cond2[1]/5),UvR1_Ylim, paste("MIGs with 1 read:", min_UvR1_Cond2, "%"))
dev.off()
###UvR2_Cond2
xlimUvR2_Cond2 <- rev(range(UvR2_Cond2$seq))
xlimUvR2_Cond2[1] <- xlimUvR2_Cond2[1]-((xlimUvR2_Cond2[1])*0.03)
xlimUvR2_Cond2[2] <- xlimUvR2_Cond2[2]+((xlimUvR2_Cond2[1])*0.03)
svg(paste(args[2],"_UMIvsReads_2read",".svg", sep = ""), width=7, height=7)  
plot(UvR2_Cond2$seq,UvR2_Cond2$count, log = "y", type="l", xlim = xlimUvR2_Cond2, ylim=c(2,UvR2_Ylim), col="dark blue", lwd=7, xlab="Ranked MIG", ylab="Log 10 (Number of Reads)", main=args[2])
text((xlimUvR2_Cond2[1]/5),UvR2_Ylim, paste("MIGs with 2 reads:", min_UvR2_Cond2, "%"))
dev.off()