#!/usr/bin/env Rscript
##############################################################################
##                           Description                                   ##
##############################################################################
#Title: Spectratype
#Description: 
#Installing and loading required packages
if(!require(ggplot2))
{install.packages("ggplot2")
}
if(!require(reshape))
{install.packages("reshape")
}
if(!require(RColorBrewer))
{install.packages("RColorBrewer")
}
library(ggplot2)
library(reshape)
library(RColorBrewer)
library(matrixStats)
library(dplyr)
##############################################################################
##                           Initial Setup                                 ##
##############################################################################
#set working directories
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#####################################

args = commandArgs(trailingOnly=TRUE)


args =c("Alpha", "2.NOD", "vdjtools.2.NOD-clones_2.txt", "spectratype_2reads_$PREFIX_NAME.svg")





#READ files in the directory:
Files_data <- as.data.frame(args[3])
colnames(Files_data) <- "files"
Files_data$files <- as.character(Files_data$files)

#Do the short name and sort
tmp_1 <- as.data.frame(gsub("-clones_2_cutted.txt", "", Files_data[,1]))
tmp_1 <- as.data.frame(gsub("vdjtools.", "", tmp_1[,1]))
colnames(tmp_1) <- c("cond")
tmp_1 <- tmp_1 %>% separate(cond, c("A", "B"), sep = "([\\.])")
Files_data$shortname <- paste(tmp_1$B, tmp_1$A, sep = "_")
Files_data$strain <- tmp_1$B
Files_data$number <- as.numeric(tmp_1$A)
Files_data <- Files_data[order(Files_data$strain, Files_data$number, decreasing = FALSE),]
Files_data$number <- NULL





###############################################################################
##                       Process Data                                        ##
###############################################################################
Data <- matrix(data=NA, nrow=8, ncol = 2)
Data[,1] <- c(10,11,12,13,14,15,16,17)
colnames(Data) <- c("Len", args[2])
  tmp_1 <- read.table(args[3], header = TRUE)
  tmp_2 <- tmp_1$count
  tmp_3 <- nchar(as.character(tmp_1$cdr3aa))
  tmp_4 <- cbind(tmp_2, tmp_3)
  colnames(tmp_4) <- c("counts", "Len")
  tmp_5 <- aggregate(counts ~ Len, data=tmp_4, sum)
  tmp_5$counts <- (tmp_5$counts/sum(tmp_5$counts))*100
  indx <- match(c(10,11,12,13,14,15,16,17), tmp_5[,1] , nomatch = 0) 
  tmp_5 <- tmp_5[indx,]
  Data[,2] <- tmp_5$counts
  
##########################################
#Process data to graph
##########################################

tmp_1 <- Data[,c(1,8,9,12)]
tmp_1 <- as.data.frame(tmp_1)
tmp_1$strain <- rep(args[2], nrow(tmp_1))
colnames(tmp_1) <- c("Len", "mean", "SD", "pValue", "strain")
tmp_2 <- Data[,c(1,10,11,12)]
tmp_2 <- as.data.frame(tmp_2)
tmp_2$strain <- rep(cond2, nrow(tmp_2))
colnames(tmp_2) <- c("Len", "mean", "SD", "pValue", "strain")
# join the samples
Data_graph <- rbind(tmp_1, tmp_2)
#change X values
Data_graph[,1] <- c(seq(from=1, to=24, by=3), seq(from=1, to=24, by=3)+1)

######### Add pValue
for (i in 1:nrow(Data_graph)) {
  if( as.numeric(Data[i,8]) < as.numeric(Data[i,10]))  {
    pvalue <- wilcox.test(as.numeric(Data[i,2:4]), as.numeric(Data[i,5:7]),  alternative = "less", paired = FALSE)
    Data[i,12] <- round(pvalue$p.value, 2) } else{
      pvalue <- wilcox.test(as.numeric(Data[i,5:7]), as.numeric(Data[i,2:4]),  alternative = "greater", paired = FALSE)
      Data[i,12] <- round(pvalue$p.value, 2)}}
Data_graph[,4] <- c(Data[,12], Data[,12])

### pValues stars
### *P < 0.05, **P < 0.01, ***P < 0.001, ****P < 0.0001 
one_star_Data <- Data_graph[which(Data_graph$pValue <= 0.05 & Data_graph$pValue > 0.01),]
one_star_Data <- one_star_Data[which(one_star_Data$strain == cond2),]
rownames(one_star_Data) <- NULL
two_star_Data <- Data_graph[which(Data_graph$pValue <= 0.01 & Data_graph$pValue > 0.001),]
two_star_Data <- two_star_Data[which(two_star_Data$sample == cond2),]
rownames(two_star_Data) <- NULL
three_star_Data <- Data_graph[which(Data_graph$pValue <= 0.001 & Data_graph$pValue > 0.0001),]
three_star_Data <- three_star_Data[which(three_star_Data$sample == cond2),]
rownames(three_star_Data) <- NULL
four_star_Data <- Data_graph[which(Data_graph$pValue <= 0.0001),]
four_star_Data <- four_star_Data[which(four_star_Data$sample == cond2),]
rownames(four_star_Data) <- NULL

#########################################
#      To graph general mean length 
#########################################
# Found the legths
lengths_data=NA
#
for (i in 1:nrow(Files_data)) {
   tmp_1 <- read.table(Files_data[i,1], header = TRUE)
   tmp_2 <- nchar(as.character(tmp_1$cdr3aa))
   lengths_data <- c(lengths_data, unique(tmp_2))
   lengths_data <- unique(lengths_data)
}
lengths_data <- sort(lengths_data)
lengths_data <- as.numeric(lengths_data)
Data_All_length <- matrix(data=0, ncol = 9, nrow=length(lengths_data))
colnames(Data_All_length) <- c("Len", Files_data$shortname, paste("mean_", args[2], sep = ""), paste("mean_", cond2, sep = ""))
Data_All_length[,1] <- lengths_data
# Add the data
for (i in 1:nrow(Files_data)) {
  tmp_1 <- read.table(Files_data[i,1], header = TRUE)
  tmp_2 <- tmp_1$count
  tmp_3 <- nchar(as.character(tmp_1$cdr3aa))
  tmp_4 <- cbind(tmp_2, tmp_3)
  colnames(tmp_4) <- c("counts", "Len")
  tmp_5 <- aggregate(counts ~ Len, data=tmp_4, sum)
  tmp_5$counts <- tmp_5$counts/sum(tmp_5$counts)
  indx <- match(tmp_5[,1] , lengths_data, nomatch = 0) 
  Data_All_length[indx,i+1] <- tmp_5$counts
  }
#means
Data_All_length[,8] <- rowMeans( Data_All_length[,2:4])
Data_All_length[,9] <- rowMeans( Data_All_length[,5:7])

#Add the length to graphic:
Data_All_length[,8] <- Data_All_length[,8]*Data_All_length[,1]
Data_All_length[,9] <- Data_All_length[,9]*Data_All_length[,1]

### To graph the mean
Data_graph_general <- matrix(data=NA, nrow=2, ncol=4)
colnames(Data_graph_general) <- c("strain", "mean", "SD", "pValue")
Data_graph_general[1:2,1] <- c(args[2], cond2)
#args[2]
tmp_1 <- as.data.frame(Data_All_length[,8])
tmp_1 <-  tmp_1[-(which(tmp_1 == 0)),]  
Data_graph_general[1,2] <- sum(tmp_1)   #/length(tmp_1)
Data_graph_general[1,3] <- sd(tmp_1)
#cond2
tmp_2 <- as.data.frame(Data_All_length[,9])
tmp_2 <-  tmp_2[-(which(tmp_2 == 0)),]  
Data_graph_general[2,2] <- sum(tmp_2)   #/length(tmp_2)
Data_graph_general[2,3] <- sd(tmp_2)
#pValue
if (sum(tmp_1)>sum(tmp_2)){pValue <- wilcox.test(tmp_1, tmp_2,  alternative = "greater", paired = FALSE)}else
{pValue <- wilcox.test(tmp_1, tmp_2,  alternative = "less", paired = FALSE)}
pValue <- pValue$p.value
Data_graph_general[1:2,4] <- c(pValue, pValue)
Data_graph_general <- as.data.frame(Data_graph_general)
#delete factors
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
Data_graph_general$strain <- as.character(Data_graph_general$strain)
Data_graph_general$mean <- as.numeric.factor (Data_graph_general$mean)
Data_graph_general$SD <- as.numeric.factor(Data_graph_general$SD)
Data_graph_general$pValue <- as.numeric.factor(Data_graph_general$pValue)
#Add positions
Data_graph_general$position <- c(1,3)
#Round the mean
Data_graph_general[,2] <-round(Data_graph_general[,2],2)

############################################################################## 
#                               Graphics                                     #
##############################################################################
##Get colors from RGB to hexadecimal:
#rgb2hex <- function(r,g,b) sprintf('#%s',paste(as.hexmode(c(r,g,b)),collapse = ''))
#rgb2hex(0,186,56) 
#c("#00ba38", "#f8766d", "#619cff")
size_chr=14
custom_theme <- function () { theme(legend.text=element_text(size=size_chr), plot.title= element_text(size=size_chr, hjust = 0.5), 
                                    axis.text.x = element_text(size=size_chr, colour = "black"), axis.text.y = element_text(size=size_chr, colour = "black"), 
                                    axis.title.x= element_text(size=size_chr+2, hjust = 0.5) , axis.title.y= element_text(size=size_chr+2),
                                    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                    panel.background = element_blank(), axis.ticks.length=unit(.25, "cm"), axis.ticks = element_line(size = 1.5),
                                    axis.line = element_line(size = 1, colour = "black"))}  
custom_theme_2 <- function () { theme(legend.text=element_text(size=size_chr), plot.title= element_text(size=size_chr, hjust = 0.5), 
                                    axis.text.x = element_text(size=size_chr, colour = "black"), axis.text.y = element_text(size=size_chr, colour = "black"), 
                                    axis.title.x= element_text(size=size_chr+2, hjust = 0.5) , axis.title.y= element_text(size=size_chr+2),
                                    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none",
                                    panel.background = element_blank(), axis.ticks.length=unit(.25, "cm"), axis.ticks.x=element_blank(),axis.ticks = element_line(size = 1.5),
                                    axis.line = element_line(size = 1, colour = "black"))}  


ggplot(Data_graph, aes(x = Len, y = mean, fill = strain))+ 
  geom_bar(width = 0.92, stat = "identity") +
  xlab("CDR3 length, aa") +
  ylab("Frequency") +
  labs(fill="", title = "") +
  scale_y_continuous(expand = c(0, 0)) +
  custom_theme() +
  guides(fill = guide_legend(reverse = FALSE))+
  expand_limits(y=c(0, 0.34))+ 
  scale_fill_manual(values=c("#156cff", "#ff492e"))+
  scale_x_continuous(breaks=seq(from=1.5, to=22.5, by=3), labels = seq(from=10, to=17))+
  geom_errorbar(aes(ymin=mean-SD, ymax=mean+SD), width=0.3, position=position_dodge(.9), size=0.7, alpha=0.6)+
  #stars
  annotate(geom="text", x=one_star_Data$Len, y=one_star_Data$mean+one_star_Data$SD+0.01, label="*",color="black", size=6)+
  annotate(geom="text", x=two_star_Data$Len, y=two_star_Data$mean+two_star_Data$SD+0.01, label="*",color="black", size=6)+
  annotate(geom="text", x=three_star_Data$Len, y=three_star_Data$mean+three_star_Data$SD+0.01, label="*",color="black", size=6)+
  annotate(geom="text", x=four_star_Data$Len, y=four_star_Data$mean+four_star_Data$SD+0.01, label="*",color="black", size=6)+
ggsave(paste("Spectratype_", args[2], "_", cond2, "_", args[1],".png", sep=""), width=14, height=10, units = "cm")



ggplot(Data_graph_general, aes(x = position, y = mean, fill = strain))+ 
  geom_bar(width = 0.8, stat = "identity") +
  xlab("") +
  ylab("CDR3 length, aa") +
  custom_theme_2() +
  expand_limits(y=c(0, 19), x=c(0,4))+ 
  scale_fill_manual(values=c("#156cff", "#ff492e"))+
  scale_y_continuous(breaks=seq(from=1, to=19, by=2), expand = c(0, 0))+
  scale_x_continuous(breaks=c(1,3), label = c("B6", "NOD"))+
  geom_errorbar(aes(ymin=mean-SD, ymax=mean+SD), width=0.3, position=position_dodge(.9), size=0.7, alpha=0.6)+
  annotate("segment", x = 1, xend = 3, y = 17, yend = 17, colour = "grey40")+
  annotate(geom="text", x=2, y=17.6, label="ns",color="grey40", size=5)+
  geom_text(aes(label=mean), vjust=4, color="white",position = position_dodge(0.9), size=3.5)+
  ggsave(paste("Spectratype_mean_", args[2], "_", cond2, "_", args[1],".png", sep=""), width=10, height=10, units = "cm")


