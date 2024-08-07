# Load the packages

library(tidyverse)
library(ggrepel)
library(ggplot2)
library(ggpubr)
library(scales)

load(here::here("ALL_eWAT_new.Rdata"))

#add a column based on another column values

ALL_eWAT_new$FDR_log <- -log10(ALL_eWAT_new$FDR)

# add a column of NAs
ALL_eWAT_new$diffexpressed <- NA

# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP"
ALL_eWAT_new$diffexpressed[ALL_eWAT_new$logFC > 0.5 & ALL_eWAT_new$FDR < 0.05] <- "Up"

# if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
ALL_eWAT_new$diffexpressed[ALL_eWAT_new$logFC < -0.5 & ALL_eWAT_new$FDR < 0.05] <- "Down"


TOP10logFCup <- c('Hephl1','Gm11295','Elovl6','1700113H08Rik','Gm50337','Fgf13',
                  'Slc2a5','Cidea','Lctl','BC049762')
TOP10logFCdown <- c('Cd7','Art4','Epha3','Col6a5','Mmrn1','Camp','Hpca','Cldn11',
                    'Ighv1-9','Crisp1')
TOP10geneslogFC <- c(TOP10logFCup, TOP10logFCdown)


TOP20logFCup <- c('Hephl1','Gm11295','Elovl6','1700113H08Rik','Gm50337','Fgf13',
                  'Slc2a5','Cidea','Lctl','BC049762','Mal2','Ntsr2','Plekhg6',
                  'Sycp3','Arhgdig','Ciart','Stpg1','Ppp1r3b','Bmp3','Gm9860')
TOP20logFCdown <- c('Rgs1','Rab15','Pirt','C7','Fndc5','Gm20619','Duox2','Prss34',
                    'Ccl17','Reln','Cd7','Art4','Epha3','Col6a5','Mmrn1','Camp','Hpca','Cldn11',
                    'Ighv1-9','Crisp1')

# Now write down the name of genes beside the points...
# Create a new column "delabel", that will contain the name of differentially expressed genes (NA in case they are not)
ALL_eWAT_new$delabel1<- NA

#find positions genes pathway and save gene names as hits
hits1 <- which(ALL_eWAT_new$Gene %in% TOP10geneslogFC)

#add names to new column
ALL_eWAT_new$delabel1[hits1] <- ALL_eWAT_new$Gene[hits1]

ALL_eWAT_new$TOP20fcup <- NA
TOP20fcup_hits <- which(ALL_eWAT_new$Gene %in% TOP20logFCup)
ALL_eWAT_new$TOP20fcup[TOP20fcup_hits] <- ALL_eWAT_new$Gene[TOP20fcup_hits]

ALL_eWAT_new$TOP20fcdown <- NA
TOP20fcdown_hits <- which(ALL_eWAT_new$Gene %in% TOP20logFCdown)
ALL_eWAT_new$TOP20fcdown[TOP20fcdown_hits] <- ALL_eWAT_new$Gene[TOP20fcdown_hits]

#volcano plot

RNAseq_eWAT_plot1 <- ALL_eWAT_new %>%
    ggplot(aes(x=logFC, y=FDR_log, colour=diffexpressed, label=delabel1)) +
    geom_point(size = 3, alpha = 0.6) +
    geom_text_repel(size = 5, nudge_y = 0.5, fontface=3) +
    labs(colour= "Differentially expressed",
         y = "-log[10](FDR)",
         x = "logFC") +
    xlim(-7,4) +
    ylim(0,6) +
    theme_pubr(base_size = 16, legend = c("top"))+
    geom_vline(xintercept=c(-0.5, 0.5), col="black", linetype=2) + # Add lines
    geom_hline(yintercept=-log10(0.05), col="black", linetype=2) +
    scale_color_manual(values=c("#6FD952","#A659D9"),breaks = c('Down', 'Up'))+
    guides(colour = guide_legend(override.aes = aes(label = "")))

RNAseq_eWAT_plot1

#remove NAs
ALL_eWAT_new_TOP20up <- ALL_eWAT_new[!is.na(ALL_eWAT_new$TOP20fcup),]

TOP20up_plot_circles_v <- ALL_eWAT_new_TOP20up %>%
    ggplot(aes(y=reorder(Gene, logFC), x = 1, colour= logFC, size=FDR_log)) +
    geom_point() +
    labs(title = "Top 20 upregulated genes",
         x = NULL,
         y = NULL,
         size = expression("-log"[10]*"(FDR)")) +
    theme_classic2(base_size = 16)+
    theme(axis.line.y = element_blank(),
          axis.text.y = element_text(hjust=1, face=3),
          axis.ticks.y = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5))+
    scale_colour_gradient2(low = muted("#6FD952"),
                           mid = "white",
                           high = muted("#A659D9"),
                           midpoint = 0)+
    guides(colour = guide_colorbar(order = 2), size = guide_legend(order = 1))
TOP20up_plot_circles_v

ALL_eWAT_new_TOP20down <- ALL_eWAT_new[!is.na(ALL_eWAT_new$TOP20fcdown),]

TOP20down_plot_circles_v <- ALL_eWAT_new_TOP20down %>%
    ggplot(aes(y=reorder(Gene, logFC), x = 1, colour= logFC, size=FDR_log)) +
    geom_point() +
    labs(title = "Top 20 downregulated genes",
         y = NULL,
         x = NULL,
         size = expression("-log"[10]*"(FDR)")) +
    theme_classic2(base_size = 16)+
    theme(axis.line.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(hjust=1, face=3),
          plot.title = element_text(hjust = 0.3))+
    scale_colour_gradient2(low = muted("#6FD952"),
                           mid = "white",
                           high = muted("#A659D9"),
                           midpoint = 0)+
    guides(colour = guide_colorbar(order = 2), size = guide_legend(order = 1))

TOP20down_plot_circles_v


