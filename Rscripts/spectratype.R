
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
#################################
#Take arguments from command line:
args = commandArgs(trailingOnly=TRUE)

args =c(Alpha, spectratype_2reads_$PREFIX_NAME.spectraV.wt.txt, test.svg)


##############################################################################
##                           Initial Setup                                 ##
##############################################################################
#set working directories
#setwd(working_directory)
#READ file:
df_B6 <- read.table(args[2], sep="\t", comment="", header = TRUE)
df_NOD <- read.table(args[3], sep="\t", comment="", header = TRUE)
#Palette selection
if (args[1]=="Alpha"){
palette <- read.table("Palette_120colors_VAsegments.txt", sep ="\t", header = TRUE, comment="", quote="")
}else{palette <- read.table("Palette_36colors_VBsegments.txt", sep ="\t", header = TRUE, comment="", quote="")}
###############################################################################
##                       Process Data and graph                              ##
###############################################################################

##############
#B6
##############
df_B6[, 1:ncol(df_B6)] <- apply(df_B6[, 1:ncol(df_B6)], 2, as.numeric)
df.m_B6 <- melt(df_B6, id = "Len")
# palette  
TRV_B6 <- as.matrix(colnames (df_B6))
TRV_B6 <- as.data.frame(TRV_B6[3:14,])
colnames(TRV_B6) <- c("TRV")
TRV_B6$TRV <- chartr(".", "-", TRV_B6$TRV)
TRV_B6$colors <- NA
m=0
for (i in TRV_B6$TRV) {
  m=m+1
  number_match <- as.numeric(which(i == palette[1])[[1]])
  color_tmp <-as.data.frame(palette[number_match,3])
  color_tmp <- as.matrix(color_tmp)
  TRV_B6 [m,2] <- color_tmp[1,1]
}
rm(number_match, color_tmp, m,i)
palette_B6 <- as.character(TRV_B6$colors)
####################################### plotting
lab_names_B6 <- chartr(".", "-", df.m_B6$variable) 
lab_names_B6 <- unique(lab_names_B6)
svg(args[4])
ggplot(df.m_B6, aes(x = Len, y = value, fill = variable)) +
  geom_bar(width = 1, stat = "identity") +
  xlab("CDR3 length, bp") +
  labs(fill="Variable segment", title = "B6") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values=c("grey75", palette_B6), label=lab_names_B6) +
  theme_bw() +
  theme(legend.text=element_text(size=8), axis.title.y=element_blank()) +
  guides(fill = guide_legend(reverse = TRUE))+
  expand_limits(y=c(0, 0.29))
dev.off()

##############
#NOD
##############
df_NOD[, 1:ncol(df_NOD)] <- apply(df_NOD[, 1:ncol(df_NOD)], 2, as.numeric)
df.m_NOD <- melt(df_NOD, id = "Len")
# palette  
TRV_NOD <- as.matrix(colnames (df_NOD))
TRV_NOD <- as.data.frame(TRV_NOD[3:14,])
colnames(TRV_NOD) <- c("TRV")
TRV_NOD$TRV <- chartr(".", "-", TRV_NOD$TRV)
TRV_NOD$colors <- NA
m=0
for (i in TRV_NOD$TRV) {
  m=m+1
  number_match <- as.numeric(which(i == palette[1])[[1]])
  color_tmp <-as.data.frame(palette[number_match,3])
  color_tmp <- as.matrix(color_tmp)
  TRV_NOD [m,2] <- color_tmp[1,1]
}
rm(number_match, color_tmp, m,i)
palette_NOD <- as.character(TRV_NOD$colors)
####################################### plotting
lab_names_NOD <- chartr(".", "-", df.m_NOD$variable) 
lab_names_NOD <- unique(lab_names_NOD)
svg(args[5])
ggplot(df.m_NOD, aes(x = Len, y = value, fill = variable)) +
  geom_bar(width = 1, stat = "identity") +
  xlab("CDR3 length, bp") +
  labs(fill="Variable segment", title = "NOD") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values=c("grey75", palette_NOD), label=lab_names_NOD) +
  theme_bw() +
  theme(legend.text=element_text(size=8), axis.title.y=element_blank()) +
  guides(fill = guide_legend(reverse = TRUE))+
  expand_limits(y=c(0, 0.29))
dev.off()


##########################################
#Process data to graph both strains
##########################################
# Remove lower lengths
to_keep <- as.data.frame(c(30,33,36,39,42,45,48,51))
colnames(to_keep) <- "numbers"
# B6
indx <- match(df.m_B6$Len, to_keep$numbers, nomatch = 0) 
indx <- as.data.frame(indx)
a <- as.numeric(which(indx$indx > 0))
df.m_B6_graph <- df.m_B6[a,]
rm(indx, a)
# NOD
indx <- match(df.m_NOD$Len, to_keep$numbers, nomatch = 0) 
indx <- as.data.frame(indx)
a <- as.numeric(which(indx$indx > 0))
df.m_NOD_graph <- df.m_NOD[a,]
#Add +1.1 to data NOD
df.m_NOD_graph_plus1 <- df.m_NOD_graph
df.m_NOD_graph_plus1$Len <- df.m_NOD_graph_plus1$Len+1.1
#add simbol to NON Vsegments:
df.m_NOD_graph_plus1_text <- data.frame(matrix(ncol = 1, nrow = 104))
colnames(df.m_NOD_graph_plus1_text) <- "variable"
for (i in df.m_NOD_graph_plus1$variable){
  number_match <- as.numeric(which(i == df.m_NOD_graph_plus1$variable))
  p <- paste(i, " ", sep="")
  for (j in number_match){
    df.m_NOD_graph_plus1_text[j, 1] <- p
  }
}
rm(p, i, j)
df.m_NOD_graph_plus1$variable  <- df.m_NOD_graph_plus1_text$variable
#join both
df_both <- rbind(df.m_B6_graph, df.m_NOD_graph_plus1)
palette_both <- c("grey75", palette_B6, "grey75", palette_NOD)

#####################################################  By strain
to_keep_plus <- as.data.frame(to_keep$numbers +1.1)
colnames(to_keep_plus) <- "numbers"

#B6
B6_Len <- matrix(data= NA, nrow=8, ncol=2)
colnames(B6_Len) <- c("Len", "value")
m=1
for (i in to_keep$numbers) {
  number_Len <- as.data.frame(df_both[which(i == df_both$Len), 3])
  colnames(number_Len) <- "Len"
  number_value_mean <- sum(number_Len$Len)
  B6_Len[m,1] <- i
  B6_Len[m,2] <- number_value_mean
  m=m+1
}
rm(m, i, number_Len, number_value_mean)

#NOD
NOD_Len <- matrix(data= NA, nrow=8, ncol=2)
colnames(NOD_Len) <- c("Len", "value")
m=1
for (i in to_keep_plus$numbers) {
  number_Len <- as.data.frame(df_both[which(i == df_both$Len), 3])
  colnames(number_Len) <- "Len"
  number_value_mean <- sum(number_Len$Len)
  NOD_Len[m,1] <- i
  NOD_Len[m,2] <- number_value_mean
  m=m+1
}
rm(m, i, number_Len, number_value_mean)

single_histogram <- as.data.frame(rbind(B6_Len, NOD_Len))
single_histogram$variable <- c(rep("B6", 8) , rep("NOD", 8))

####################################### plotting
lab_names_1 <- chartr(".", "-", df.m_B6_graph$variable) 
lab_names_1 <- unique(lab_names_1)
lab_names_2 <- chartr(".", "-", df.m_NOD_graph_plus1$variable) 
lab_names_2 <- unique(lab_names_2)
lab_names_both <- c(lab_names_1, lab_names_2)

# B6 and NOD
svg(args[6] , width=12.5, height=9)
ggplot(df_both, aes(x = Len, y = value, fill = variable)) +
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values=palette_both, label=lab_names_both) +
  xlab(" B6   +                   +                  +                +                 +                  +                  +                  + \n NOD        +                 +                 +                  +                  +                 +                  +                  + \n                                                                              CDR3 length, bp") +
  ylab("Frequency") +
  labs(fill="       Variable segment \n\n B6                   NOD ", title = "B6 and NOD") +
  scale_y_continuous(expand = c(0, 0)) +  #, label=palette_both$TRV) +
  theme_bw() +
  theme(legend.text=element_text(size=8), plot.title= element_text(size=20, hjust = 0.5), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14), axis.title.x= element_text(size=14, hjust = 0) , axis.title.y= element_text(size=14)) +
  #guides(fill = guide_legend(reverse = TRUE))+
  expand_limits(y=c(0, 0.29))+ 
  scale_x_continuous(breaks=c(c(30.55,33.55,36.55,39.55,42.55,45.55,48.55,51.55)), labels = c("30","33","36","39","42","45","48","51"))
dev.off()

### By strain
svg(paste(args[6], "_single_color.svg", sep = "") , width=12.5, height=9)
ggplot(single_histogram, aes(x = Len, y = value, fill = variable))+ 
  geom_bar(width = 1, stat = "identity") +
  xlab("CDR3 length, bp") +
  ylab("Frequency") +
  labs(fill="", title = "B6 and NOD") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(legend.text=element_text(size=14), plot.title= element_text(size=20, hjust = 0.5), axis.text.x = element_text(size=14), axis.text.y = element_text(size=14), axis.title.x= element_text(size=14, hjust = 0.5) , axis.title.y= element_text(size=14)) +
  guides(fill = guide_legend(reverse = FALSE))+
  expand_limits(y=c(0, 0.29))+ 
  scale_x_continuous(breaks=c(c(30.55,33.55,36.55,39.55,42.55,45.55,48.55,51.55)), labels = c("30","33","36","39","42","45","48","51"))+
  scale_fill_manual(values=c("blue1", "brown1"))
dev.off()