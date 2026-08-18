#Initialisation Projet
####Open package ####
if (!require("pacman")) install.packages("pacman")

pacman::p_load(
  tidyverse, here, patchwork, ggpubr, reshape2, paletteer, 
  RColorBrewer, gridExtra, cowplot,
  magick, grid, rio, stringr, pheatmap, ggnewscale, rstatix,renv,ggnewscale

)

#### Style ####
palette_TDG <- c(
  "Fresh" = "#0072B2",
  "Thaw" = "#009E73",
  "Purified" = "#CC79A7",
  "Non Purified" = "#0072B2"
)

palette_TDG_pastel <- c(
  "Fresh" = "#A8CBE3",
  "Thaw" = "#97C1A9",
  "Purified" = "#C9A9D3"
)
palette_ICP <- c(
  "BTLA" = "#3B5B92",
  "CD16" = "#B45F5F",
  "CD39" = "#4F8A7A",
  "LAG-3" = "#8A6BAE",
  "PD-1" = "#C58A3A",
  "TIGIT" = "#5E7F9F",
  "TIM-3" = "#A65E86"
)

palette_diff <- c(
  "Naive \n CD45RA+ CCR7+" = "#3B5B92",
  "Central memory \n CD45RA- CCR7+" = "#4F8A7A",
  "Effector memory \n CD45RA- CCR7-" = "#C58A3A",
  "Terminal effector memory \n CD45RA+ CCR7-" = "#A65E86"
)

theme_blood <- function() {
  theme_classic() +  # Fond classique et épuré
    theme(
      axis.text.x = element_text(size = 6, face = "bold", color = "black"),  # Labels des axes X
      axis.text.y = element_text(size = 6, face = "bold", color = "black"),
      axis.title = element_text(size = 8, face = "bold"),
      axis.title.y = element_text(size = 8, face = "bold", color = "black"),  # Labels des axes X
      legend.title = element_text(size = 7, face = "bold"),  # Titre de la légende
      legend.text = element_text(size = 6,face="bold"),
      legend.box.spacing = unit(0.1, "cm"),
      legend.key.size = unit(0.4, "cm"),  
      plot.margin = margin(t = 0.6,  # Top margin
                           r = 0.5,  # Right margin
                           b = 0,  # Bottom margin
                           l = 1), # Left margin 
      plot.title = element_text(hjust = 0.5, size = 9, face = "bold"),  # Titre du graphique centré
      strip.placement = "outside",  # Strips à l'extérieur
      strip.background = element_blank(),  # Suppression du cadre autour des strips
      strip.text = element_text(size = 6, face = "bold", angle = 0)
    )
}