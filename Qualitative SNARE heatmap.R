

#import necessary packages
library(gtools)
library(plyr)
library(readxl)
library(dplyr)
library(ggplot2)
library(reshape2)
library(scales)
library(grid)

#read source files
source_data <- read_excel("src/qualitative input.xlsx")
source_data <- source_data[mixedorder(source_data$Name),]
#clear empty rows
source_data_clean <- source_data[rowSums(is.na(source_data[3:12])) != ncol(source_data[3:12]),]
#remove empty columns
source_data_clean <- source_data_clean[colSums(is.na(source_data_clean)) != nrow(source_data_clean)]

#Change type column to factor for colouring by row
source_data_clean$type <- as.factor(source_data_clean$type)

#replace NA by 0
source_data_clean[is.na(source_data_clean)] <- 0

#Fix column names
source_data_clean$Name <- tolower(source_data_clean$Name)
source_data_clean$Name <- paste0(toupper(substr(source_data_clean$Name, 1, 1)), substr(source_data_clean$Name, 2, nchar(source_data_clean$Name)))

#transform data to facilitate heatmap processing
source_data_clean.m <- melt(source_data_clean)
source_data_clean.s <- ddply(source_data_clean.m, .(variable), transform, rescale = scale(value))

#reverse y axis labels (were z-a, now a-z)
flevels <- levels(as.factor(source_data_clean.s$Name))
flevels <- mixedsort(flevels, decreasing = TRUE)

#define colours per SNARE type
source_data_clean.s$valueoffset <- source_data_clean.s$value + 100*(as.numeric(source_data_clean.s$type)-1)
scalerange <- range(source_data_clean.s$value)
gradientends <- scalerange+ rep(c(0,100,200,300,400), each = 2)
colorends <- c("white", "firebrick2", "white", "orange2", "white", "forestgreen", "white", "mediumaquamarine", "white", "royalblue2")

pdf("Qualitative phagosome proteomics heatmap.pdf", paper = "a4")
ggplot(source_data_clean.s, aes(variable, Name)) +
  geom_tile(aes(fill = valueoffset, shape = substr(type, 1, 3)), colour = "lightgrey", show.legend = TRUE) +
  geom_vline(xintercept = 2.5) +
  scale_fill_gradientn(colours = colorends, values = rescale(gradientends), guide = FALSE) +
  scale_x_discrete("", expand = c(0, 0)) +
  scale_y_discrete("", limits = flevels, expand = c(0, 0)) +
  theme_grey(base_size = 14) +
  theme(legend.position = "right",
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0, face = "bold")) +
  guides(shape = guide_legend("SNARE type", keyheight = 2, override.aes = list(fill = c("firebrick2", "orange2", "forestgreen", "mediumaquamarine", "royalblue2"))))
dev.off()
