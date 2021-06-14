#!/usr/bin/env Rscript
##############################################################################
##                           Description                                   ##
##############################################################################
#Title: Circoplot
#Description: 
#Installing and loading required packages
if(!require(circlize))
{install.packages("circlize")
}
if(!require(RColorBrewer))
{install.packages("RColorBrewer")
}
if(!require(stringr))
{install.packages("stringr")
}
library(circlize)
library(RColorBrewer)
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
df_B6 <- read.table(args[2], sep="\t", comment="")
df_NOD <- read.table(args[3], sep="\t", comment="")

#Palette selection
if (args[1]=="Alpha"){
  palette_V <- read.table("Palette_120colors_VAsegments.txt", sep ="\t", header = TRUE, comment="", quote="")
  palette_J <- read.table("Palette_61colors_JAsegments.txt", sep ="\t", header = TRUE, comment="", quote="")
}else{palette_V <- read.table("Palette_36colors_VBsegments.txt", sep ="\t", header = TRUE, comment="", quote="")
      palette_J <- read.table("Palette_14colors_JBsegments.txt", sep ="\t", header = TRUE, comment="", quote="")}
###############################################################################
##                          Process Data ALL SAMPLES                         ##
###############################################################################
#B6
n <- nrow(df_B6)
m <- ncol(df_B6)
rn = as.character(df_B6[2:n,1])
cn = apply(df_B6[1,2:m], 2 , as.character)
mat_B6 <- matrix(apply(df_B6[2:n, 2:m], 1:2, as.numeric), n - 1, m-1) * 100
n <- nrow(df_B6)
m <- ncol(df_B6)
duplicates <- intersect(rn, cn)
rownames(mat_B6) <- replace(rn, rn==duplicates, paste("V", duplicates, sep=""))
colnames(mat_B6) <- replace(cn, cn==duplicates, paste("J", duplicates, sep=""))
rm(n, m, rn, cn)
#NOD
n <- nrow(df_NOD)
m <- ncol(df_NOD)
rn = as.character(df_NOD[2:n,1])
cn = apply(df_NOD[1,2:m], 2 , as.character)
mat_NOD <- matrix(apply(df_NOD[2:n, 2:m], 1:2, as.numeric), n - 1, m-1) * 100
n <- nrow(df_NOD)
m <- ncol(df_NOD)
duplicates <- intersect(rn, cn)
rownames(mat_NOD) <- replace(rn, rn==duplicates, paste("V", duplicates, sep=""))
colnames(mat_NOD) <- replace(cn, cn==duplicates, paste("J", duplicates, sep=""))
rm(n, m, rn, cn)

####  sort
#B6
col_sum = apply(mat_B6, 2, sum)
row_sum = apply(mat_B6, 1, sum)
mat_B6 <- mat_B6[order(row_sum), order(col_sum)]
#NOD
col_sum = apply(mat_NOD, 2, sum)
row_sum = apply(mat_NOD, 1, sum)
mat_NOD <- mat_NOD[order(row_sum), order(col_sum)]

#### equal number of characters for visualizaiton
#B6
rn_B6 <- rownames(mat_B6)
cn_B6 <- colnames(mat_B6)
maxrn <- max(nchar(rn_B6))
maxcn <- max(nchar(cn_B6))
for(i in seq_len(length(rn_B6))) {
  rn_B6[i] <- paste(rn_B6[i], paste(rep(" ", maxrn - nchar(rn_B6[i])), collapse = ''))
}
for(i in seq_len(length(cn_B6))) {
  cn_B6[i] <- paste(cn_B6[i], paste(rep(" ", maxcn - nchar(cn_B6[i])), collapse = ''))
}
rm(maxrn, maxcn)

#NOD
rn_NOD <- rownames(mat_NOD)
cn_NOD <- colnames(mat_NOD)
maxrn <- max(nchar(rn_NOD))
maxcn <- max(nchar(cn_NOD))
for(i in seq_len(length(rn_NOD))) {
  rn_NOD[i] <- paste(rn_NOD[i], paste(rep(" ", maxrn - nchar(rn_NOD[i])), collapse = ''))
}
for(i in seq_len(length(cn_NOD))) {
  cn_NOD[i] <- paste(cn_NOD[i], paste(rep(" ", maxcn - nchar(cn_NOD[i])), collapse = ''))
}
rm(maxrn, maxcn)

####Delete Spaces in the rownames and colnames and replace dash simbol
#B6
rn_B6 <- str_replace_all(rownames(mat_B6), "\\s", "")
rn_B6 <- str_replace_all(rn_B6, "/", "-")
cn_B6 <- str_replace_all(colnames(mat_B6), "\\s", "")
cn_B6 <- str_replace_all(cn_B6, "/", "-")
rownames(mat_B6) <- rn_B6
colnames(mat_B6) <- cn_B6
#NOD
rn_NOD <- str_replace_all(rownames(mat_NOD), "\\s", "")
rn_NOD <- str_replace_all(rn_NOD, "/", "-")
cn_NOD <- str_replace_all(colnames(mat_NOD), "\\s", "")
cn_NOD <- str_replace_all(cn_NOD, "/", "-")
rownames(mat_NOD) <- rn_NOD
colnames(mat_NOD) <- cn_NOD
####Re-order the table:
mat_B6_tmp <- as.data.frame(mat_B6)
mat_NOD_tmp <- as.data.frame(mat_NOD)

####Vsegments
#B6
mat_B6_tmp <- t(mat_B6_tmp)
mat_B6_tmp <- as.data.frame(mat_B6_tmp)
mat_B6_tmp$orther_row <- NA
m=0 
V_tmp=0
for (i in row.names(mat_B6_tmp)) {
  m=m+1
  number_match <- which(i == palette_V[1])[[1]]
  V_tmp <-as.numeric (palette_V[number_match,4])
  #  V_tmp <- as.matrix(V_tmp)
  mat_B6_tmp [m,ncol(mat_B6_tmp)] <- V_tmp
}
rm(number_match, V_tmp, m,i)
mat_B6_tmp <- mat_B6_tmp[ order(mat_B6_tmp$orther_row), ]
mat_B6_tmp$orther_row <- NULL

#NOD
mat_NOD_tmp <- t(mat_NOD_tmp)
mat_NOD_tmp <- as.data.frame(mat_NOD_tmp)
mat_NOD_tmp$orther_row <- NA
m=0 
V_tmp=0
for (i in row.names(mat_NOD_tmp)) {
  m=m+1
  number_match <- which(i == palette_V[1])[[1]]
  V_tmp <-as.numeric (palette_V[number_match,4])
  #  V_tmp <- as.matrix(V_tmp)
  mat_NOD_tmp [m,ncol(mat_NOD_tmp)] <- V_tmp
}
rm(number_match, V_tmp, m,i)
mat_NOD_tmp <- mat_NOD_tmp[ order(mat_NOD_tmp$orther_row), ]
mat_NOD_tmp$orther_row <- NULL

####J segments
#B6
mat_B6_tmp <- t(mat_B6_tmp)
mat_B6_tmp <- as.data.frame(mat_B6_tmp)
mat_B6_tmp$orther_row <- NA
m=0 
V_tmp=0
for (i in row.names(mat_B6_tmp)) {
  m=m+1
  number_match <- which(i == palette_J[1])[[1]]
  V_tmp <-as.numeric (palette_V[number_match,4])
  #  V_tmp <- as.matrix(V_tmp)
  mat_B6_tmp [m,ncol(mat_B6_tmp)] <- V_tmp
}
rm(number_match, V_tmp, m,i)
mat_B6_tmp <- mat_B6_tmp[ order(mat_B6_tmp$orther_row), ]
mat_B6_tmp$orther_row <- NULL
mat_B6 <- as.matrix(mat_B6_tmp)

#NOD
mat_NOD_tmp <- t(mat_NOD_tmp)
mat_NOD_tmp <- as.data.frame(mat_NOD_tmp)
mat_NOD_tmp$orther_row <- NA
m=0 
V_tmp=0
for (i in row.names(mat_NOD_tmp)) {
  m=m+1
  number_match <- which(i == palette_J[1])[[1]]
  V_tmp <-as.numeric (palette_V[number_match,4])
  #  V_tmp <- as.matrix(V_tmp)
  mat_NOD_tmp [m,ncol(mat_NOD_tmp)] <- V_tmp
}
rm(number_match, V_tmp, m,i)
mat_NOD_tmp <- mat_NOD_tmp[ order(mat_NOD_tmp$orther_row), ]
mat_NOD_tmp$orther_row <- NULL
##############################################################################
##                          Graphics                                        ##
##############################################################################

#######Palette  
#### J segments
#B6
jnames_B6 <- as.data.frame(rownames(mat_B6))
jnames_B6$colors <- NA
m=0
for (i in jnames_B6$`rownames(mat_B6)`) {
  m=m+1
  number_match <- as.numeric(which(i == palette_J[1])[[1]])
  color_tmp <-as.data.frame(palette_J[number_match,3])
  color_tmp <- as.matrix(color_tmp)
  jnames_B6 [m,2] <- color_tmp[1,1]
}
rm(number_match, color_tmp, m,i)
#NOD
jnames_NOD <- as.data.frame(rownames(mat_NOD))
jnames_NOD$colors <- NA
m=0
for (i in jnames_NOD$`rownames(mat_NOD)`) {
  m=m+1
  number_match <- as.numeric(which(i == palette_J[1])[[1]])
  color_tmp <-as.data.frame(palette_J[number_match,3])
  color_tmp <- as.matrix(color_tmp)
  jnames_NOD [m,2] <- color_tmp[1,1]
}
rm(number_match, color_tmp, m,i)

#### V segments
#B6
vnames_B6 <- as.data.frame(colnames(mat_B6))
vnames_B6$colors <- NA
m=0
for (i in vnames_B6$`colnames(mat_B6)`) {
  m=m+1
  number_match <- as.numeric(which(i == palette_V[1])[[1]])
  color_tmp <-as.data.frame(palette_V[number_match,3])
  color_tmp <- as.matrix(color_tmp)
  vnames_B6 [m,2] <- color_tmp[1,1]
}
rm(number_match, color_tmp, m,i)
#NOD
vnames_NOD <- as.data.frame(colnames(mat_NOD))
vnames_NOD$colors <- NA
m=0
for (i in vnames_NOD$`colnames(mat_NOD)`) {
  m=m+1
  number_match <- as.numeric(which(i == palette_V[1])[[1]])
  color_tmp <-as.data.frame(palette_V[number_match,3])
  color_tmp <- as.matrix(color_tmp)
  vnames_NOD [m,2] <- color_tmp[1,1]
}
rm(number_match, color_tmp, m,i)

#### Names
#B6
rcols_B6 <- as.character(jnames_B6$colors)
names(rcols_B6) <- rownames(mat_B6)
ccols_B6 <- as.character(vnames_B6$colors)
names(ccols_B6) <- colnames(mat_B6)
#NOD
rcols_NOD <- as.character(jnames_NOD$colors)
names(rcols_NOD) <- rownames(mat_NOD)
ccols_NOD <- as.character(vnames_NOD$colors)
names(ccols_NOD) <- colnames(mat_NOD)

#### Plot B6
svg(args[4])
circos.par(gap.degree = c(rep(1, nrow(mat_B6)-1), 10, rep(1, ncol(mat_B6)-1), 15), start.degree = 5)
chordDiagram(mat_B6, annotationTrack = "grid",
             grid.col = c(rcols_B6, ccols_B6), title("B6"),
             preAllocateTracks = list(track.height = 0.2), transparency = 0.5)
circos.trackPlotRegion(track.index = 1, bg.border = NA,
                       panel.fun = function(x, y) {
                         sector.name = get.cell.meta.data("sector.index")
                         xlim = get.cell.meta.data("xlim")
                         ylim = get.cell.meta.data("ylim")
                         circos.text(mean(xlim), ylim[1], cex = 0.5, sector.name, facing = "clockwise", adj = c(0, 0.5))
                       }
)

circos.clear()
dev.off()

#### Plot NOD
svg(args[5])
circos.par(gap.degree = c(rep(1, nrow(mat_NOD)-1), 10, rep(1, ncol(mat_NOD)-1), 15), start.degree = 5)
chordDiagram(mat_NOD, annotationTrack = "grid",
             grid.col = c(rcols_NOD, ccols_NOD), title("NOD"),
             preAllocateTracks = list(track.height = 0.2), transparency = 0.5)
circos.trackPlotRegion(track.index = 1, bg.border = NA,
                       panel.fun = function(x, y) {
                         sector.name = get.cell.meta.data("sector.index")
                         xlim = get.cell.meta.data("xlim")
                         ylim = get.cell.meta.data("ylim")
                         circos.text(mean(xlim), ylim[1], cex = 0.5, sector.name, facing = "clockwise", adj = c(0, 0.5))
                       }
)

circos.clear()
dev.off()