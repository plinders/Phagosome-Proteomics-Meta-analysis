



#import necessary packages
library(gtools)
library(plyr)
library(readxl)
library(dplyr)
library(ggplot2)
library(reshape2)
library(scales)
library(grid)
library(data.table)
library(extrafont)
library(gridExtra)
library(gtable)

font_import(pattern = "[C/c]alibri")

#read source files
source_data <-
  read_excel("src/Final meta-analysis scaled by row_new.xlsx")
source_data <- source_data[mixedorder(source_data$`Gene names`),]
source_data[which(source_data$`Gene names` == "Stx5;Stx5a"), 1] <-
  "Stx5"
snare_list_source <- read_excel("src/human SNAREs.xlsx")
snare_list <- as.character(snare_list_source$Name)


#filter source data for only required SNAREs
hmap_input <-
  dplyr::filter(source_data,
                tolower(source_data$`Gene names`) %in% tolower(snare_list))
source_data$`Gene names` <- tolower(source_data$`Gene names`)
snare_list_new <- dplyr::select(snare_list_source, Name, type)
snare_list_new$Name <- tolower(snare_list_new$Name)
new_hmap_input <-
  merge(snare_list_new, source_data, by.x = "Name", by.y = "Gene names")

write.csv(new_hmap_input, "hmap_input.csv")

#Fix column names
new_hmap_input$Name <-
  paste0(toupper(substr(new_hmap_input$Name, 1, 1)), substr(new_hmap_input$Name, 2, nchar(new_hmap_input$Name)))

new_hmap_input <- rename(new_hmap_input, GosR1 = Gosr1)

new_hmap_input$Name[which(new_hmap_input$Name == "Gosr1")] <-
  "GosR1"
new_hmap_input$Name[which(new_hmap_input$Name == "Stxbp5")] <-
  "STXBP5"
new_hmap_input$Name[which(new_hmap_input$Name == "Snap23")] <-
  "SNAP23"
new_hmap_input$Name[which(new_hmap_input$Name == "Snap29")] <-
  "SNAP29"
new_hmap_input$Name[which(new_hmap_input$Name == "Vamp3")] <-
  "VAMP3"
new_hmap_input$Name[which(new_hmap_input$Name == "Vamp4")] <-
  "VAMP4"
new_hmap_input$Name[which(new_hmap_input$Name == "Vamp7")] <-
  "VAMP7"
new_hmap_input$Name[which(new_hmap_input$Name == "Vamp8")] <-
  "VAMP8"




# diff_df <- cbind(Name = new_hmap_input$Name, as.data.frame(new_hmap_input[14:24] - new_hmap_input[3:13]))
# diff_df_mean <- diff_df %>% mutate(mean = rowMeans(diff_df[2:12], na.rm = TRUE))
# diff_df_mean <- diff_df_mean %>% select(Name, mean)
# 
# diff_df_mean.m <- melt(diff_df_mean)


#Change type column to factor for colouring by row
new_hmap_input$type <- as.factor(new_hmap_input$type)
new_hmap_input$type = factor(new_hmap_input$type, levels(new_hmap_input$type)[c(1, 2, 4, 3, 5)])
#remove unnecessary columns from Excel analyses
new_hmap_input <-
  dplyr::select(new_hmap_input, Name, type, starts_with("20"))
#sort columns alphabetically
new_hmap_input <-
  new_hmap_input[, c(names(new_hmap_input)[1:2], sort(names(new_hmap_input)[3:13]), sort(names(new_hmap_input)[14:26]))]
#reorder columns by type
new_hmap_input <- arrange(new_hmap_input, type)
new_hmap_input <- filter(new_hmap_input, Name != "STXBP5")

new_hmap_input <- new_hmap_input %>% mutate(spacer = 0) %>% mutate(Difference = 0)

#transform data to facilitate heatmap processing
new_hmap_input.m <- melt(new_hmap_input)
#mean_diff.m <- new_hmap_input.m %>% mutate(type = replace(as.character(type), variable == "mean", "S_Mean")) %>% as.data.frame()
new_hmap_input.s <-
  ddply(new_hmap_input.m, .(variable), transform, rescale = scale(value))

#reverse y axis labels (were z-a, now a-z)
flevelsdf <- new_hmap_input.s[1:2]

flevels_qa <- dplyr::filter(flevelsdf, type == "Qa")
flevels_qb <- dplyr::filter(flevelsdf, type == "Qb")
flevels_qc <- dplyr::filter(flevelsdf, type == "Qc")
flevels_qbc <- dplyr::filter(flevelsdf, type == "Qbc")
flevels_r <- dplyr::filter(flevelsdf, type == "R")

flevels_qa <- levels(as.factor(flevels_qa$Name))
flevels_qb <- levels(as.factor(flevels_qb$Name))
flevels_qc <- levels(as.factor(flevels_qc$Name))
flevels_qbc <- levels(as.factor(flevels_qbc$Name))
flevels_r <- levels(as.factor(flevels_r$Name))

flevels_qa <- mixedsort(flevels_qa, decreasing = TRUE)
flevels_qb <- mixedsort(flevels_qb, decreasing = TRUE)
flevels_qc <- mixedsort(flevels_qc, decreasing = TRUE)
flevels_qbc <- mixedsort(flevels_qbc, decreasing = TRUE)
flevels_r <- mixedsort(flevels_r, decreasing = TRUE)

flevels <-
  c(flevels_r, flevels_qbc, flevels_qc, flevels_qb, flevels_qa)



#define colours per SNARE type
new_hmap_input.s$valueoffset <-
  new_hmap_input.s$value + 100 * (as.numeric(new_hmap_input.s$type) - 1)
scalerange <- range(new_hmap_input.s$value)
gradientends <- scalerange + rep(c(0, 100, 200, 300, 400), each = 2)
colorends <-
  c("white",
    "red",
    "white",
    "green4",
    "white",
    "green",
    "white",
    "green3",
    "white",
    "blue")

pdf("170313 Quantitative phagosome proteomics heatmap.pdf", paper = "a4")

new_hmap_input.s$variable2 <-
  reorder(new_hmap_input.s$variable,
          as.numeric(new_hmap_input.s$type))

#x_tick_label <- as.data.frame(new_hmap_input.s$variable)


ggplot(new_hmap_input.s, aes(variable2, Name)) +
  geom_tile(aes(fill = valueoffset, shape = flevelsdf$type),
            colour = "lightgrey",
            show.legend = TRUE) +
  geom_vline(xintercept = 11.5) +
  geom_hline(yintercept = 6.5, xmin = -1) +
  geom_hline(yintercept = 8.5, xmin = -1) +
  geom_hline(yintercept = 10.5, xmin = -1) +
  geom_hline(yintercept = 13.5, xmin = -1) +
  scale_fill_gradientn(colours = colorends,
                       values = rescale(gradientends),
                       guide = FALSE) +
  scale_x_discrete(
    "",
    expand = c(0, 0),
    position = "top",
    labels = c(
      "[1]",
      "[2]",
      "[3]",
      "[4]",
      "[5]",
      "[6]",
      "[7]",
      "[8]",
      "[9]",
      "[10]",
      "[11]",
      "[1]",
      "[2]",
      "[3]",
      "[4]",
      "[5]",
      "[6]",
      "[7]",
      "[8]",
      "[9]",
      "[10]",
      "[11]",
      "",
      "Change"
    )
  ) +
  scale_y_discrete("", limits = flevels, expand = c(0, 0)) +
  # theme_grey(base_size = 11) +
  theme(
    legend.position = "right",
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 0,
      vjust = 1.0
    )
  ) +
  #theme(plot.margin = unit(c(1, 1, 1, 1), "lines")) +
  annotation_custom(
    grob = textGrob(
      label = "<0.5 hr after uptake",
      hjust = 0,
      vjust = 0,
      gp = gpar(fontsize = 10)
    ),
    ymin = 25,
    xmin = -16.5
  ) +
  annotation_custom(
    grob = textGrob(
      label = ">1 hr after uptake",
      hjust = 0,
      vjust = 0,
      gp = gpar(fontsize = 10)
    ),
    ymin = 25,
    xmin = 6
  ) +
  guides(shape = guide_legend(
    "SNARE type",
    keyheight = 1,
    override.aes = list(fill = c("red", "green4", "green", "green3", "blue"))
  )) +
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    axis.text = element_text(colour = "black")
  )

# p_mean <- ggplot(mean_diff.m, aes(variable, Name)) +
#   geom_tile(aes(fill = value), colour = "lightgrey") +
#   scale_fill_gradientn(guide = guide_legend(title = "Early v late", keyheight = 1, reverse = TRUE), colours = c("blue", "white", "red"), values = rescale(range(mean_diff.m$value)), limits = c(-0.5, 0.5)) +
#  # theme(plot.margin = unit(c(1, 1, 1, 1), "lines")) +
#   scale_y_discrete("", limits = flevels, expand = c(0, 0)) +
#   scale_x_discrete(
#     "",
#     expand = c(0, 0),
#     position = "top",
#     labels = c(
#       "Difference"
#     )
#   ) +
#   theme(
#     legend.position = "right",
#     axis.ticks = element_blank(),
#     axis.text.x = element_text(
#       angle = 90,
#       hjust = 0,
#       vjust = 1.0
#     ),
#     axis.text.y = element_blank()
#   )
# grid.arrange(p, p_mean, ncol = 2, nrow = 1, widths = c(5/6, 1/6), heights = 1)
#   
# gA <- ggplotGrob(p)
# gB <- ggplotGrob(p_mean)
# 
# grid.arrange(gA, gB, widths = c(5/6, 1/6))
# 
# gA$heights <- gB$heights
# grid.draw(gA)
# grid.newpage()
# grid.draw(arrangeGrob(gA, gB, widths = c(3/4, 1/4)))

# gt <- ggplot_gtable(ggplot_build(p))
# gt$layout$clip[gt$layout$name == "panel"] <- "off"
# grid.draw(gt)
# 
# gt_mean <- ggplot_gtable(ggplot_build(p_mean))
# gt$layout$clip[gt$layout$name == "panel"] <- "off"
# grid.draw(gt_mean)

dev.off()

ggsave(filename = "final.png",
       plot = grid.draw(arrangeGrob(gB, gA, widths = c(1/6, 5/6))),
       dpi = 300)

ggsave(filename = "170214 Quantitative phagosome proteomics heatmap.png",
      plot = gt,
     dpi = 300)

ggsave(filename = "170214 Quantitative phagosome proteomics heatmap_mean.png",
       plot = gt_mean,
       dpi = 300)

