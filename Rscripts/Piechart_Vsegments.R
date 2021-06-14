#!/usr/bin/env Rscript
##############################################################################
##                           Description                                   ##
##############################################################################
#Title: Circoplot
#Description: 
#Installing and loading required packages
#if(!require(grDevices))
#{install.packages("grDevices")
#}
if(!require(plyr))
{install.packages("plyr")
}
if(!require(stringr))
{install.packages("stringr")
}
#library(grDevices)
library(plyr)
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
df_B6 <- read.table(args[2], sep="\t", comment="", header = TRUE)
df_NOD <- read.table(args[3], sep="\t", comment="", header = TRUE)

#Palette selection
if (args[1]=="Alpha"){
  palette_V <- read.table("Palette_120colors_VAsegments.txt", sep ="\t", header = TRUE, comment="", quote="")
  palette_J <- read.table("Palette_61colors_JAsegments.txt", sep ="\t", header = TRUE, comment="", quote="")
}else{palette_V <- read.table("Palette_36colors_VBsegments.txt", sep ="\t", header = TRUE, comment="", quote="")
palette_J <- read.table("Palette_14colors_JBsegments.txt", sep ="\t", header = TRUE, comment="", quote="")}
###############################################################################
##                          Process Data ALL SAMPLES                         ##
###############################################################################
####Take V segments and his counts
df_B6 <- df_B6[,c(1,5)]
df_NOD <- df_NOD[,c(1,5)]

####Change symbols in TCR segments:
#B6
df_B6$v <- str_replace_all(df_B6$v, "/", "-")
#NOD
df_NOD$v <- str_replace_all(df_NOD$v, "/", "-")
#Consolidate data by sum
df_B6 <- ddply(df_B6,"v",numcolwise(sum))
df_NOD <- ddply(df_NOD,"v",numcolwise(sum))

####Add missing segments
#B6
B6a <- as.character(palette_V$V_segments) 
B6b <- as.character(df_B6$v)
differents_B6 <- as.data.frame(setdiff(B6a,B6b))
differents_B6$count <- seq(from = 0, to = 0, length.out = nrow(differents_B6))
colnames(differents_B6) <- c("v", "count")
df_B6 <- rbind(df_B6, differents_B6)
rm(B6a, B6b, differents_B6)
#NOD
NODa <- as.character(palette_V$V_segments)
NODb <- as.character(df_NOD$v)
differents_NOD <- as.data.frame(setdiff(NODa,NODb))
differents_NOD$count <- seq(from = 0, to = 0, length.out = nrow(differents_NOD))
colnames(differents_NOD) <- c("v", "count")
df_NOD <- rbind(df_NOD, differents_NOD)
rm(NODa, NODb, differents_NOD)

####Re-order the table:
#B6
df_B6$orther_row <- NA
m=0 
V_tmp=0
for (i in df_B6$v) {
  m=m+1
  number_match <- which(i == palette_V[1])[[1]]
  V_tmp <-as.numeric (palette_V[number_match,4])
  df_B6 [m,ncol(df_B6)] <- V_tmp
}
rm(number_match, V_tmp, m,i)
df_B6 <- df_B6[ order(df_B6$orther_row), ]
df_B6$orther_row <- NULL
#NOD
df_NOD$orther_row <- NA
m=0 
V_tmp=0
for (i in df_NOD$v) {
  m=m+1
  number_match <- which(i == palette_V[1])[[1]]
  V_tmp <-as.numeric (palette_V[number_match,4])
  df_NOD [m,ncol(df_NOD)] <- V_tmp
}
rm(number_match, V_tmp, m,i)
df_NOD <- df_NOD[ order(df_NOD$orther_row), ]
df_NOD$orther_row <- NULL
# Add 1 to all to avoid missing V segments
df_B6$count <- df_B6$count +1
df_NOD$count <- df_NOD$count +1

#Add Short name column
df_B6$short <- gsub("TRAV", "", df_B6$v, fixed = TRUE)
df_NOD$short <- gsub("TRAV", "", df_NOD$v, fixed = TRUE)

###############################################################################
##                          Function pie  chart                              ##
###############################################################################

#graphics::pie
pie_modif <-function (x, labels = names(x), edges = 200, radius = 0.8, clockwise = FALSE, 
                      init.angle = if (clockwise) 90 else 0, density = NULL, angle = 45, 
                      col = NULL, border = NULL, lty = NULL, main = NULL, ...) 
{
  if (!is.numeric(x) || any(is.na(x) | x < 0)) 
    stop("'x' values must be positive.")
  if (is.null(labels)) 
    labels <- as.character(seq_along(x))
  else labels <- as.graphicsAnnot(labels)
  x <- c(0, cumsum(x)/sum(x))
  dx <- diff(x)
  nx <- length(dx)
  plot.new()
  pin <- par("pin")
  xlim <- ylim <- c(-1, 1)
  if (pin[1L] > pin[2L]) 
    xlim <- (pin[1L]/pin[2L]) * xlim
  else ylim <- (pin[2L]/pin[1L]) * ylim
  dev.hold()
  on.exit(dev.flush())
  plot.window(xlim, ylim, "", asp = 1)
  if (is.null(col)) 
    col <- if (is.null(density)) 
      c("white", "lightblue", "mistyrose", "lightcyan", 
        "lavender", "cornsilk")
  else par("fg")
  if (!is.null(col)) 
    col <- rep_len(col, nx)
  if (!is.null(border)) 
    border <- rep_len(border, nx)
  if (!is.null(lty)) 
    lty <- rep_len(lty, nx)
  angle <- rep(angle, nx)
  if (!is.null(density)) 
    density <- rep_len(density, nx)
  twopi <- if (clockwise) 
    -2 * pi
  else 2 * pi
  t2xy <- function(t) {
    t2p <- twopi * t + init.angle * pi/180
    list(x = radius * cos(t2p), y = radius * sin(t2p), an=t2p)
  }
  for (i in 1L:nx) {
    n <- max(2, floor(edges * dx[i]))
    P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
    polygon(c(P$x, 0), c(P$y, 0), density = density[i], angle = angle[i], 
            border = border[i], col = col[i], lty = lty[i])
    P <- t2xy(mean(x[i + 0:1]))
    lab <- as.character(labels[i])
    if (!is.na(lab) && nzchar(lab)) {
      lines(c(1, 1.05) * P$x, c(1, 1.05) * P$y)
      text(1.1 * P$x, 1.1 * P$y, labels[i], xpd = TRUE, 
           srt = ifelse(P$x < 0, P$an/pi*180+180, P$an/pi*180),
           adj = ifelse(P$x < 0, 1, 0), ...)
    }
  }
  title(main = main, ...)
  invisible(NULL)
}

###############################################################################
##                                   Graphics                                ##
###############################################################################
#Colors
palette <-as.character(palette_V$colors)  #distinctColorPalette(100)
#text size
par(ps=9)

#B6
svg(args[4], width=9, height=9)
pie_modif(df_B6$count, labels = df_B6$short, border = NA, main= paste("B6"," V segments") , col=palette,
    init.angle = 90, clockwise = TRUE) 
dev.off()

#NOD
svg(args[5],  width=9, height=9)
pie_modif(df_NOD$count, labels = df_NOD$short, border = NA, main= paste("NOD"," V segments") , col=palette,
          init.angle = 90, clockwise = TRUE) 
dev.off()