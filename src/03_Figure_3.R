####Figure 3 Papier Léa TGD ####
source("./src/00_init_project.R")

#### Comparaison Depletion ####
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
         UTplusT_normalized_depletion_by_UTminusT = as.numeric(UTplusT_normalized_depletion_by_UTminusT))


Df_depletion_UT <- Df_depletion |> 
  filter(UTplusT_normalized_depletion_by_UTminusT > 0) |> 
  select(Donor,Thaw,Purification,UTplusT_normalized_depletion_by_UTminusT,Purif)

##### Depletion Thaw #####
Df_depletion_UT_thaw <- Df_depletion_UT |> 
  filter(Purification != "Purified") |> 
  mutate(Donor = str_remove(Donor, "\\*$")) 

plot_depletion_UT_thaw <- ggplot(Df_depletion_UT_thaw, 
                                 aes(x = Purification, y = UTplusT_normalized_depletion_by_UTminusT)) +
  geom_boxplot(aes(fill = Purification),width = 0.6, outliers = FALSE,
               linewidth = 0.7,
               color = "black") +
  scale_fill_manual(values = palette_TDG) +
  ggnewscale::new_scale_color() +
  geom_jitter(aes(color = Purification), width = 0.20, height = 0,size = 1.7, alpha = 0.75) +
  scale_color_manual(values = palette_TDG_pastel) +
  stat_compare_means(method = "wilcox.test",aes(group = Purification), hide.ns = FALSE,
                     paired = FALSE,label = "p.signif", size = 3, vjust = 1, hjust = 0.5, label.x.npc = "middle") +
  theme_blood() +
  theme(strip.text = element_text(size = 8, face = "bold"),
        legend.position =  "none") +
  labs(x = "", y = "Normalized depletion \n by UT- T cells")

##### Depletion Purif #####
plot_depletion_UT_Purif <- ggplot(Df_depletion_UT, aes(x = Purif, y = UTplusT_normalized_depletion_by_UTminusT)) +
  geom_boxplot(aes(fill = Purif),width = 0.6, outliers = FALSE,
               linewidth = 0.7,
               color = "black") +
  scale_fill_manual(values = palette_TDG) +
  ggnewscale::new_scale_color() +
  geom_jitter(aes(color = Purification), width = 0.20, height = 0,size = 1.7, alpha = 0.75) +
  scale_color_manual(values = palette_TDG_pastel) +
  stat_compare_means(method = "wilcox.test",aes(group = Purif), hide.ns = FALSE,
                     paired = FALSE,label = "p.signif", size = 3, vjust = 1, hjust = 0.5, label.x.npc = "middle") +
  theme_blood() +
  theme(strip.text = element_text(size = 8, face = "bold"),
        legend.position =  "none") +
  labs(x = "", y = "Normalized depletion \n by UT- T cells")

##### Depletion Purif Paired #####
Df_depletion_UT_paired <- Df_depletion_UT |> 
  group_by(Donor) |>
  filter(n_distinct(Purif) == n_distinct(Df_depletion_UT$Purif)) |>
  ungroup()




