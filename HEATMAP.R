# Load the packages
source(here::here("C:/Users/lmq835/Desktop/All CBMR stuff/1 InFood my project/InFood2/RNA seq eWAT/eWAT InFood2.0/Heatmap paper/Heatmap_eWAT/package-loading.R"))

heat <- read_excel("C:/Users/lmq835/Desktop/All CBMR stuff/1 InFood my project/InFood2/RNA seq eWAT/eWAT InFood2.0/Heatmap paper/Heatmap_eWAT/HEATMAP log_CPM_corrected_gene_name_by_sample.xlsx",
                   sheet = "HEATMAP")

# change column name for x column print(df)
colnames(heat)[1] <- "Genes"

#remove specific columns

heat <- heat[ -c(22) ]

#select rows with DEGs
## open excel with DEGs list
deg <- read_excel("C:/Users/lmq835/Desktop/All CBMR stuff/1 InFood my project/InFood2/RNA seq eWAT/eWAT InFood2.0/Heatmap paper/Heatmap_eWAT/HEATMAP log_CPM_corrected_gene_name_by_sample.xlsx",
                  sheet = "DEG")

#filter rows with DEGs

heat_DEGs <- heat %>%
  filter(Genes %in% deg$DEG)

#order rows to keep DEG order

orderDEG <- c(deg$DEG)

heat_DEGs <- heat_DEGs[na.omit(match(orderDEG, heat_DEGs$Genes)), ]

#To assign row names to the genes names
rownames(heat_DEGs) <- heat_DEGs$Genes

#subset columns by group

columns <- c("0159_1", "0159_2", "0159_3", "0159_4", "0159_5", "0159_11",
             "0159_12", "0159_13", "0159_14", "0159_15", "0159_6", "0159_7", "0159_8",
             "0159_9", "0159_10", "0159_16", "0159_17", "0159_18", "0159_19", "0159_20")

#Convert to data matrix, then remove column with genes

DM_heat_DEGs <- data.matrix(heat_DEGs, rownames.force = heat_DEGs$Genes)
DM_heat_DEGs <- DM_heat_DEGs[, columns]

#annotate names of groups

# Define groups for columns
groups <- c("Control", "Control", "Control", "Control", "Control", "Control",
            "Control", "Control", "Control", "Control", "Food Insecurity", "Food Insecurity", "Food Insecurity", "Food Insecurity",
            "Food Insecurity", "Food Insecurity", "Food Insecurity", "Food Insecurity", "Food Insecurity",
            "Food Insecurity")

# Assign colours to the groups
group_colors <- c("Control" = "#BF8C4E", "Food Insecurity" = "#087F8C")
col_colors <- group_colors[groups]

#Functions needed to generate heatmap https://www.datanovia.com/en/lessons/heatmap-in-r-static-and-interactive-visualization/
library("gplots")

#Heatmap with some changes in margins, and no clustering by columns (samples)

my_colors <- colorRampPalette(c("#6FD952", "white", "#A659D9"))(100)

heatmap.2(DM_heat_DEGs, scale = "row", col = my_colors, ColSideColors = col_colors,
          trace = "none", key.title = "", density.info = "none",
          lhei = c(1,5), lwid = NULL, dendrogram = "none", Colv=FALSE, Rowv=FALSE, labCol = "")


while (!is.null(dev.list())) dev.off()

