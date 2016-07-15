

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
source_data <- read_excel("src/160704 Proteomics meta-analysis scaled by row.xlsx")
source_data <- source_data[mixedorder(source_data$`Gene names`),]
source_data[which(source_data$`Gene names` == "Stx5;Stx5a"), 1] <- "Stx5"
snare_list_source <- read_excel("src/human SNAREs.xlsx")
snare_list <- as.character(snare_list_source$Name)


#filter source data for only required SNAREs
hmap_input <- dplyr::filter(source_data, tolower(source_data$`Gene names`) %in% tolower(snare_list))
source_data$`Gene names` <- tolower(source_data$`Gene names`)
snare_list_new <- dplyr::select(snare_list_source, Name, type)
snare_list_new$Name <- tolower(snare_list_new$Name)
new_hmap_input <- merge(snare_list_new, source_data, by.x = "Name", by.y = "Gene names")
#Fix column names
new_hmap_input$Name <- paste0(toupper(substr(new_hmap_input$Name, 1, 1)), substr(new_hmap_input$Name, 2, nchar(new_hmap_input$Name)))

#Change type column to factor for colouring by row  
new_hmap_input$type <- as.factor(new_hmap_input$type)
#remove unnecessary columns from Excel analyses
new_hmap_input <- dplyr::select(new_hmap_input, Name, type, starts_with("20"))
#sort columns alphabetically
new_hmap_input <- new_hmap_input[,c(names(new_hmap_input)[1:2], sort(names(new_hmap_input)[3:12]), sort(names(new_hmap_input)[13:23]))]

#transform data to facilitate heatmap processing
new_hmap_input.m <- melt(new_hmap_input)
new_hmap_input.s <- ddply(new_hmap_input.m, .(variable), transform, rescale = scale(value))

#reverse y axis labels (were z-a, now a-z)
flevels <- levels(as.factor(new_hmap_input.s$Name))
flevels <- mixedsort(flevels, decreasing = TRUE)

#define colours per SNARE type
new_hmap_input.s$valueoffset <- new_hmap_input.s$value + 100*(as.numeric(new_hmap_input.s$type)-1)
scalerange <- range(new_hmap_input.s$value)
gradientends <- scalerange+ rep(c(0,100,200,300,400), each = 2)
colorends <- c("white", "firebrick2", "white", "orange2", "white", "forestgreen", "white", "mediumaquamarine", "white", "royalblue2")

pdf("Quantitative phagosome proteomics heatmap.pdf", paper = "a4")


ggplot(new_hmap_input.s, aes(variable, Name)) +
  geom_tile(aes(fill = valueoffset, shape = substr(type, 1, 3)), colour = "lightgrey", show.legend = TRUE) +
  geom_vline(xintercept = 10.5) +
  scale_fill_gradientn(colours = colorends, values = rescale(gradientends), guide = FALSE) +
  scale_x_discrete("", expand = c(0, 0)) +
  scale_y_discrete("", limits = flevels, expand = c(0, 0)) +
  theme_grey(base_size = 11) +
  theme(legend.position = "right",
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0, face = "bold")) +
  guides(shape = guide_legend("SNARE type", keyheight = 2, override.aes = list(fill = c("firebrick2", "orange2", "forestgreen", "mediumaquamarine", "royalblue2"))))

dev.off()