set.seed(54645)

####Figure 3 Papier Léa TGD ####
source("./src/00_init_project.R")

#### Depletion GA ####
Df_depletion <- import("./data/DF/260512_Normalized_depletion_All.xlsx") |> 
  mutate(Donor = recode(Donor,"BC240112A*_1" = "BC240112A_1*",
                        "BC230331A*_1" = "BC230331A_1*",
                        "BC231020A*_1" = "BC231020A_1*",
                        "BC240322A*_1" = "BC240322A_1*",
                        "BC240209B*_1" = "BC240209B_1*")) |> 
  mutate(Thaw = if_else(str_detect(Donor, "\\*$"), "Thaw", "Fresh")) |> 
  mutate(Purif = if_else(Purified == "Yes", "Purified", "Non Purified"))   |> 
  mutate(Purification = if_else(Purif == "Purified", "Purified", Thaw),
         Purification = factor(Purification, levels = c("Fresh","Thaw","Purified")),
         UTplusT_normalized_depletion_by_UTminusT = as.numeric(UTplusT_normalized_depletion_by_UTminusT),
         GAminusT_normalized_depletion_by_UTminusT = as.numeric(GAminusT_normalized_depletion_by_UTminusT),
         GAplusT_normalized_depletion_by_UTminusT = as.numeric(GAplusT_normalized_depletion_by_UTminusT)) |> 
  rename( "RL + Tγδ UT" = "UTplusT_normalized_depletion_by_UTminusT", 
          "RL - Tγδ GA101" = "GAminusT_normalized_depletion_by_UTminusT",
          "RL + Tγδ GA101"="GAplusT_normalized_depletion_by_UTminusT")

##### Non Purif vs Purif ####
Df_depletion_longer <- Df_depletion |> 
  pivot_longer(cols = c("RL + Tγδ UT","RL - Tγδ GA101","RL + Tγδ GA101"),
               names_to = "Conditions",
               values_to = "Per_Depletion") |> 
  mutate(Conditions = factor(Conditions, levels =  c("RL + Tγδ UT","RL - Tγδ GA101","RL + Tγδ GA101")))


plot_depletion_GA101 <- ggplot(Df_depletion_longer, aes(x = Purif, y = Per_Depletion)) +
  facet_wrap(~Conditions, strip.position = "bottom")+
  stat_summary(aes(color = Purif),fun = mean,geom = "crossbar", width = 0.6,
               linewidth = 1.2,
               fatten = 1) +
  scale_color_manual(values = palette_TDG) +
  ggnewscale::new_scale_color() +
  geom_beeswarm(aes(color = Purification), cex = 4, corral = "wrap", corral.width = 0.6, size = 1.7, alpha = 0.75) +
  scale_color_manual(values = palette_TDG_pastel) +
  stat_compare_means(method = "wilcox.test",aes(group = Purif), hide.ns = FALSE,
                     paired = FALSE,label = "p.signif", size = 3, vjust = 1, hjust = 0.5, label.x.npc = "middle") +
  theme_blood() +
  theme(strip.text = element_text(size = 8, face = "bold"),
        legend.position =  "none") +
  labs(x = "", y = "Normalized depletion \n by UT- T cells")


##### Plot All samples ####
test_comp <- list(
  c("RL + Tγδ UT", "RL - Tγδ GA101"),
  c("RL + Tγδ UT", "RL + Tγδ GA101"),
  c("RL - Tγδ GA101", "RL + Tγδ GA101")
)

plot_depletion_GA101_all <- ggplot(Df_depletion_longer, aes(x = Conditions, y = Per_Depletion)) +
  stat_summary(color = "#B45F5F",fun = mean,geom = "crossbar", width = 0.6,
               linewidth = 1.2,
               fatten = 1) +
  scale_color_manual(values = palette_TDG) +
  ggnewscale::new_scale_color() +
  geom_beeswarm(aes(color = Purification), cex = 4, corral = "wrap", corral.width = 0.6, size = 1.7, alpha = 0.75) +
  scale_color_manual(values = palette_TDG_pastel) +
  stat_compare_means(method = "wilcox.test", comparisons = test_comp, p.adjust.method = "BH", hide.ns = FALSE,
                     paired = TRUE,label = "p.signif", size = 3, vjust = 0, hjust = 0.5, label.x.npc = "middle") +
  theme_blood() +
  theme(strip.text = element_text(size = 8, face = "bold"),
        legend.position =  "none") +
  labs(x = "", y = "Normalized depletion \n by UT- T cells")

##### Non Purif Only ####
Df_depletion_GA101_non_purif <- Df_depletion_longer |> 
  filter(Purif == "Non Purified")


test_comp <- list(
  c("RL + Tγδ UT", "RL - Tγδ GA101"),
  c("RL + Tγδ UT", "RL + Tγδ GA101"),
  c("RL - Tγδ GA101", "RL + Tγδ GA101")
)

plot_depletion_GA101_non_purif <- ggplot(Df_depletion_GA101_non_purif, aes(x = Conditions, y = Per_Depletion)) +
  stat_summary(aes(color = Purif),fun = mean,geom = "crossbar", width = 0.6,
               linewidth = 1.2,
               fatten = 1) +
  scale_color_manual(values = palette_TDG) +
  ggnewscale::new_scale_color() +
  geom_beeswarm(aes(color = Purification), cex = 2, corral = "wrap", corral.width = 0.6, size = 1.7, alpha = 0.75) +
  scale_color_manual(values = palette_TDG_pastel) +
  stat_compare_means(method = "wilcox.test", comparisons = test_comp, p.adjust.method = "BH", hide.ns = FALSE,
                     paired = TRUE,label = "p.signif", size = 3, vjust = 0, hjust = 0.5, label.x.npc = "middle") +
  theme_blood() +
  theme(strip.text = element_text(size = 8, face = "bold"),
        legend.position =  "none") +
  labs(x = "", y = "Normalized depletion \n by UT- T cells")

##### Purif Only ####
Df_depletion_GA101_purif <- Df_depletion_longer |> 
  filter(Purif == "Purified")


test_comp <- list(
  c("RL + Tγδ UT", "RL - Tγδ GA101"),
  c("RL + Tγδ UT", "RL + Tγδ GA101"),
  c("RL - Tγδ GA101", "RL + Tγδ GA101")
)

plot_depletion_GA101_purif <- ggplot(Df_depletion_GA101_purif, aes(x = Conditions, y = Per_Depletion)) +
  stat_summary(aes(color = Purif),fun = mean,geom = "crossbar", width = 0.6,
               linewidth = 1.2,
               fatten = 1) +
  scale_color_manual(values = palette_TDG) +
  ggnewscale::new_scale_color() +
  geom_beeswarm(aes(color = Purification), cex = 2, corral = "wrap", corral.width = 0.6, size = 1.7, alpha = 0.75) +
  scale_color_manual(values = palette_TDG_pastel) +
  stat_compare_means(method = "wilcox.test", comparisons = test_comp, p.adjust.method = "BH", hide.ns = FALSE,
                     paired = TRUE,label = "p.signif", size = 3, vjust = 0, hjust = 0.5, label.x.npc = "middle") +
  theme_blood() +
  theme(strip.text = element_text(size = 8, face = "bold"),
        legend.position =  "none") +
  labs(x = "", y = "Normalized depletion \n by UT- T cells")

##### FC Block ####
Liste_Fc_block <- Df_depletion |> 
  filter(FcBlock_reel == "Yes") |> 
  group_by(Donor) |> 
  distinct(Donor) |> 
  pull(Donor)

Df_depletion_longer_Fc_Block <- Df_depletion |> 
  filter(Donor %in% Liste_Fc_block) |> 
  mutate(FcBlock_reel = recode(FcBlock_reel, "No" = "No FC Block", "Yes" = "FC Block")) |> 
  mutate(FcBlock_reel = factor(FcBlock_reel, levels = c("No FC Block","FC Block"))) |> 
  mutate(GAplusT_normalized_depletion_by_UTplusT = as.numeric(GAplusT_normalized_depletion_by_UTplusT),
         GAplusT_normalized_depletion_by_GAminusT = as.numeric(GAplusT_normalized_depletion_by_GAminusT)) |> 
  pivot_longer(cols = c("RL + Tγδ UT","RL + Tγδ GA101"),
               names_to = "Conditions",
               values_to = "Per_Depletion") |> 
  mutate(Conditions = factor(Conditions, levels =  c("RL + Tγδ UT","RL - Tγδ GA101","RL + Tγδ GA101")))
  
  
  # pivot_longer(cols = c("GAplusT_normalized_depletion_by_UTplusT","GAplusT_normalized_depletion_by_GAminusT"),
  #              names_to = "Conditions",
  #              values_to = "Per_Depletion")


plot_depletion_GA101_Fcblock <- ggplot(Df_depletion_longer_Fc_Block, aes(x = FcBlock_reel, y = Per_Depletion)) +
  facet_wrap(~Conditions, strip.position = "bottom")+
  geom_line(aes(group = Donor), alpha = 0.15) +
  stat_summary(aes(color = Purif),fun = mean,geom = "crossbar", width = 0.6,
               linewidth = 1.2,
               fatten = 1) +
  scale_color_manual(values = palette_TDG) +
  ggnewscale::new_scale_color() +
  geom_point(aes(color = Purification),size = 1.7, alpha = 0.75) +
  scale_color_manual(values = palette_TDG_pastel) +
  stat_compare_means(method = "wilcox.test",aes(group = FcBlock_reel), hide.ns = FALSE,
                     paired = TRUE,label = "p.signif", size = 3, vjust = 1, hjust = 0.5, label.x.npc = "middle") +
  theme_blood() +
  theme(strip.text = element_text(size = 8, face = "bold"),
        legend.position =  "none") +
  labs(x = "", y = "Normalized depletion \n by UT- T cells")


#### CD107a ####
#Preparation
Df_cd107a <- import("./data/DF/260520_Datas_CD107a_Fcblock_Purification_tidy.xlsx") |>
  mutate(Donor = recode(Donor,"BC240112A*_1" = "BC240112A_1*",
                        "BC230331A*_1" = "BC230331A_1*",
                        "BC231020A*_1" = "BC231020A_1*",
                        "BC240322A*_1" = "BC240322A_1*",
                        "BC240209B*_1" = "BC240209B_1*")) |> 
  mutate(Thaw = if_else(str_detect(Donor, "\\*$"), "Thaw", "Fresh")) |> 
  mutate(Purif = if_else(Purified == "Yes", "Purified", "Non Purified"))   |> 
  mutate(Purification = if_else(Purif == "Purified", "Purified", Thaw),
         Purification = factor(Purification, levels = c("Fresh","Thaw","Purified"))) |> 
  mutate(across(c(UTplusT, UT_T_only, GAplusT, GA_Tonly), as.numeric)) |> 
  rename( "Tγδ UT"= "UT_T_only",
          "Tγδ GA101"= "GA_Tonly",
          "RL + Tγδ UT" = "UTplusT",
          "RL + Tγδ GA101" = "GAplusT")

##### Non Purif vs Purif #####
Df_cd107a_GA101 <- Df_cd107a |> 
  filter(!is.na(`RL + Tγδ UT`))|>
  select(-`RL + Tγδ UT`) |> 
  pivot_longer(cols = c(where(is.numeric)),
               names_to = "Conditions",
               values_to = "Per_CD107a") |> 
  mutate(Conditions = factor(Conditions, levels = c("Tγδ UT","Tγδ GA101","RL + Tγδ UT","RL + Tγδ GA101")))

#Plot
plot_CD107a_GA101 <- ggplot(Df_cd107a_GA101, aes(x = Purif, y = Per_CD107a)) +
  facet_wrap(~Conditions, strip.position = "bottom", ncol = 4)+
  stat_summary(aes(color = Purif),fun = mean,geom = "crossbar", width = 0.6,
               linewidth = 1.2,
               fatten = 1) +
  scale_color_manual(values = palette_TDG) +
  ggnewscale::new_scale_color() +
  geom_beeswarm(aes(color = Purification), cex = 4, corral = "wrap", corral.width = 0.6, size = 1.7, alpha = 0.75) +
  scale_color_manual(values = palette_TDG_pastel) +
  stat_compare_means(method = "wilcox.test",aes(group = Purif), hide.ns = FALSE,
                     paired = FALSE,label = "p.signif", size = 3, vjust = 1, hjust = 0.5, label.x.npc = "middle") +
  theme_blood() +
  theme(strip.text = element_text(size = 8, face = "bold"),
        legend.position =  "none") +
  labs(x = "", y = "% CD107a on \n CD3+ Tγδ")

##### Non Purif #####
Df_cd107a_GA101_non_purif <- Df_cd107a_GA101 |> 
  filter(Purif == "Non Purified")

test_comp <- list(
  c("Tγδ UT", "Tγδ GA101"),
  c("Tγδ UT", "RL + Tγδ GA101"),
  c("Tγδ GA101", "RL + Tγδ GA101")
)

plot_CD107a_GA101_non_purif <- ggplot(Df_cd107a_GA101_non_purif, aes(x = Conditions, y = Per_CD107a)) +
  stat_summary(aes(color = Purif),fun = mean, geom = "crossbar", width = 0.6,
               linewidth = 1.2,
               fatten = 1) +
  scale_color_manual(values = palette_TDG) +
  ggnewscale::new_scale_color() +
  geom_beeswarm(aes(color = Purification), cex = 2, corral = "wrap", corral.width = 0.6, size = 1.7, alpha = 0.75) +
  scale_color_manual(values = palette_TDG_pastel) +
  stat_compare_means(method = "wilcox.test", comparisons = test_comp, p.adjust.method = "BH", hide.ns = FALSE,
                     paired = TRUE, label = "p.signif", size = 3,hjust = 0.5, label.x.npc = "middle") +
  theme_blood() +
  theme(strip.text = element_text(size = 8, face = "bold"),
        legend.position =  "none") +
  labs(x = "", y = "% CD107a on \n CD3+ Tγδ")

##### Purif #####
Df_cd107a_GA101_purif <- Df_cd107a_GA101 |> 
  filter(Purif == "Purified")

test_comp <- list(
  c("Tγδ UT", "Tγδ GA101"),
  c("Tγδ UT", "RL + Tγδ GA101"),
  c("Tγδ GA101", "RL + Tγδ GA101")
)

plot_CD107a_GA101_purif <- ggplot(Df_cd107a_GA101_purif, aes(x = Conditions, y = Per_CD107a)) +
  stat_summary(aes(color = Purif),fun = mean, geom = "crossbar", width = 0.6,
               linewidth = 1.2,
               fatten = 1) +
  scale_color_manual(values = palette_TDG) +
  ggnewscale::new_scale_color() +
  geom_beeswarm(aes(color = Purification), cex = 2, corral = "wrap", corral.width = 0.6, size = 1.7, alpha = 0.75) +
  scale_color_manual(values = palette_TDG_pastel) +
  stat_compare_means(method = "wilcox.test", comparisons = test_comp, p.adjust.method = "BH", hide.ns = FALSE,
                     paired = TRUE, label = "p.signif", size = 3,hjust = 0.5, label.x.npc = "middle") +
  theme_blood() +
  theme(strip.text = element_text(size = 8, face = "bold"),
        legend.position =  "none") +
  labs(x = "", y = "% CD107a on \n CD3+ Tγδ")


#### Granzyme B ####
Df_granzyme_B <- import("./data/DF/260406_Datas_GrzB_Fcblock_Purification.xlsx") |>
  mutate(Donor = recode(Donor,"BC240112A*_1" = "BC240112A_1*",
                        "BC230331A*_1" = "BC230331A_1*",
                        "BC231020A*_1" = "BC231020A_1*",
                        "BC240322A*_1" = "BC240322A_1*",
                        "BC240209B*_1" = "BC240209B_1*")) |> 
  mutate(Thaw = if_else(str_detect(Donor, "\\*$"), "Thaw", "Fresh")) |> 
  mutate(Purif = if_else(Purified == "Yes", "Purified", "Non Purified"))   |> 
  mutate(Purification = if_else(Purif == "Purified", "Purified", Thaw),
         Purification = factor(Purification, levels = c("Fresh","Thaw","Purified"))) |> 
  mutate(across(c(UTminusT,UTplusT, GAplusT), as.numeric)) |> 
  rename( "RL - Tγδ UT" = "UTminusT",
          "RL + Tγδ UT" = "UTplusT",
          "RL + Tγδ GA101" = "GAplusT")

##### Non Purif vs Purif #####
Df_granzyme_B_GA101 <- Df_granzyme_B |> 
  filter(!is.na(`RL + Tγδ UT`))|>
  pivot_longer(cols = c(where(is.numeric)),
               names_to = "Conditions",
               values_to = "Conc_GrzmB") |> 
  mutate(Conditions = factor(Conditions, levels = c("RL - Tγδ UT","RL + Tγδ UT","RL + Tγδ GA101")))

#Plot
plot_Granzyme_B_GA101 <- ggplot(Df_granzyme_B_GA101, aes(x = Purif, y = Conc_GrzmB)) +
  facet_wrap(~Conditions, strip.position = "bottom", ncol = 4)+
  stat_summary(aes(color = Purif),fun = mean,geom = "crossbar", width = 0.6,
               linewidth = 1.2,
               fatten = 1) +
  scale_color_manual(values = palette_TDG) +
  ggnewscale::new_scale_color() +
  geom_beeswarm(aes(color = Purification), cex = 4, corral = "wrap", corral.width = 0.6,
                size = 1.7, alpha = 0.75) +
  scale_color_manual(values = palette_TDG_pastel) +
  stat_compare_means(method = "wilcox.test",aes(group = Purif), hide.ns = FALSE,
                     paired = FALSE,label = "p.signif", size = 3, vjust = 1, hjust = 0.5, label.x.npc = "middle") +
  theme_blood() +
  theme(strip.text = element_text(size = 8, face = "bold"),
        legend.position =  "none") +
  labs(x = "", y = "Granzyme B (pg/mL)")

##### Non Purif #####
Df_granzyme_B_GA101_non_purif <- Df_granzyme_B_GA101 |> 
  filter(Purif == "Non Purified")

test_comp <- list(
  c("RL - Tγδ UT", "RL + Tγδ UT"),
  c("RL - Tγδ UT", "RL + Tγδ GA101"),
  c("RL + Tγδ UT", "RL + Tγδ GA101")
)

plot_granzyme_B_GA101_non_purif <- ggplot(Df_granzyme_B_GA101_non_purif, aes(x = Conditions, y = Conc_GrzmB)) +
  stat_summary(aes(color = Purif),fun = mean, geom = "crossbar", width = 0.6,
               linewidth = 1.2,
               fatten = 1) +
  scale_color_manual(values = palette_TDG) +
  ggnewscale::new_scale_color() +
  geom_beeswarm(aes(color = Purification), cex = 2, corral = "wrap", corral.width = 0.6, size = 1.7, alpha = 0.75) +
  scale_color_manual(values = palette_TDG_pastel) +
  stat_compare_means(method = "wilcox.test", comparisons = test_comp, p.adjust.method = "BH", hide.ns = FALSE,
                     paired = TRUE, label = "p.signif", size = 3,hjust = 0.5, label.x.npc = "middle") +
  theme_blood() +
  theme(strip.text = element_text(size = 8, face = "bold"),
        legend.position =  "none") +
  labs(x = "", y = "Granzyme B (pg/mL)")

##### Purif #####
Df_granzyme_B_GA101_purif <- Df_granzyme_B_GA101 |> 
  filter(Purif == "Purified")

plot_granzyme_B_GA101_purif <- ggplot(Df_granzyme_B_GA101_purif, aes(x = Conditions, y = Conc_GrzmB)) +
  stat_summary(aes(color = Purif),fun = mean, geom = "crossbar", width = 0.6,
               linewidth = 1.2,
               fatten = 1) +
  scale_color_manual(values = palette_TDG) +
  ggnewscale::new_scale_color() +
  geom_beeswarm(aes(color = Purification), cex = 2, corral = "wrap", corral.width = 0.6, size = 1.7, alpha = 0.75) +
  scale_color_manual(values = palette_TDG_pastel) +
  stat_compare_means(method = "wilcox.test", comparisons = test_comp, p.adjust.method = "BH", hide.ns = FALSE,
                     paired = TRUE, label = "p.signif", size = 3,hjust = 0.5, label.x.npc = "middle") +
  theme_blood() +
  theme(strip.text = element_text(size = 8, face = "bold"),
        legend.position =  "none") +
  labs(x = "", y = "Granzyme B (pg/mL)")



#### Montage ####
##### Main ####
Figure_4 <- plot_grid(
  plot_grid(plot_depletion_GA101,plot_CD107a_GA101,legend_fig,
            ncol = 3, labels = c("A","B"), vjust = 1.1,rel_widths = c(1.1,0.9,0.2), align = "hv", axis = "bt"),
  plot_grid(plot_Granzyme_B_GA101,plot_depletion_GA101_Fcblock,
            ncol = 3, labels = c("C","D","","",""), rel_widths = c(1.1,0.6,0.5),
            vjust = 0.85), # Deuxième ligne
  nrow = 2,
  rel_heights = c(1,1),
  axis = "bt",
  align = "hv"
)

print(Figure_4)

ggsave(
  filename = "./outputs/Figures/Figure_Lea_4.png",
  plot = Figure_4,
  device = "png",
  width = 32, # largeur A4 en cm
  height = 14, # hauteur A4 en cm
  units = "cm",
  dpi = 600,
  bg = "white"
)

##### Suppl #####
Figure_4_suppl <- plot_grid(
  plot_grid(plot_depletion_GA101_non_purif,plot_depletion_GA101_purif,legend_fig,
            ncol = 3, labels = c("A","B"), vjust = 1.1,rel_widths = c(1.1,0.9,0.2), align = "hv", axis = "bt"),
  plot_grid(plot_CD107a_GA101_non_purif,plot_CD107a_GA101_purif,legend_fig,
            ncol = 3, labels = c("C","D"), vjust = 1.1,rel_widths = c(1.1,0.9,0.2), align = "hv", axis = "bt"),
  plot_grid(plot_granzyme_B_GA101_non_purif,plot_granzyme_B_GA101_purif,legend_fig,
            ncol = 3, labels = c("E","F"), vjust = 1.1,rel_widths = c(1.1,0.9,0.2), align = "hv", axis = "bt"),
  nrow = 3,
  rel_heights = c(1,1,1),
  axis = "bt",
  align = "hv"
)

print(Figure_4_suppl)

ggsave(
  filename = "./outputs/Figures/Figure_Lea_4_suppl.png",
  plot = Figure_4_suppl,
  device = "png",
  width = 34, # largeur A4 en cm
  height = 16, # hauteur A4 en cm
  units = "cm",
  dpi = 600,
  bg = "white"
)
