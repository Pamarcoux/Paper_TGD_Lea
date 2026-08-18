####Figure 2 Papier Léa TGD ####
source("./src/00_init_project.R")

#### Comparaison Thaw Fresh ####
##### % d1 or d2 ####
Df_Per_d1_d2_comparaison <- import("./data/DF/260817_Datas_CD3_d1_or_d2.ods") |> 
  mutate(Donor = recode(Donor,"BC240112A*_1" = "BC240112A_1*")) |> 
  mutate(Thaw = if_else(str_detect(Donor, "\\*$"), "Thaw", "Fresh")) |> 
  pivot_longer(cols = where(is.numeric),
               names_to = "Type_delta",
               values_to = "Per_delta") |> 
  mutate(Type_delta = recode(Type_delta, "%CD3 d1" = "CD3+ δ1+", "%CD3 d2" = "CD3+ δ2+")) |> 
  mutate(Donor = str_remove(Donor, "\\*$"))

plot_Per_d1_d2_comparison <- ggplot(Df_Per_d1_d2_comparaison, aes(x = Thaw, y = Per_delta)) +
  facet_wrap(~Type_delta, strip.position = "top" )+
  geom_jitter(aes(color = Thaw), width = 0.20, height = 0,size = 1.7, alpha = 0.75) +
  scale_color_manual(values = palette_TDG_pastel) +
  ggnewscale::new_scale_color() +
  stat_summary(aes(color = Thaw),fun = mean,geom = "crossbar", width = 0.6,
    linewidth = 1.2,
    fatten = 1) +
  scale_color_manual(values = palette_TDG) +
  stat_compare_means(method = "wilcox.test",aes(group = Thaw), hide.ns = FALSE,
    paired = FALSE,label = "p.signif", size = 3, vjust = 1, hjust = 0.5, label.x.npc = "middle") +
  theme_blood() +
  theme(strip.text = element_text(size = 8, face = "bold"),
        legend.position =  "none") +
  labs(x = "", y = "% of CD3γ9 cells at day of use",)

##### Comparison ICP #### 
#Prep
Df_icp_on_use <- import("./data/DF/260414_Datas_phenotyping_%_RFI_tidy.xlsx") |> 
  mutate(Donor = recode(Donor,"BC240112A*_1" = "BC240112A_1*",
                        "BC230331A*_1" = "BC230331A_1*",
                        "BC231020A*_1" = "BC231020A_1*",
                        "BC240322A*_1" = "BC240322A_1*",
                        "BC240209B*_1" = "BC240209B_1*")) |> 
  mutate(Thaw = if_else(str_detect(Donor, "\\*$"), "Thaw", "Fresh")) 

Df_icp_on_use_comparaison_paired <- Df_icp_on_use |>
  mutate(Donor = str_remove(Donor, "\\*$")) |> 
  filter(Purif == FALSE) |> 
  group_by(Donor) |>
  filter(n_distinct(Thaw) == n_distinct(Df_icp_on_use$Thaw)) |>
  ungroup() |> 
  mutate(across(c(CD161, CD16, CD39, `TIM-3`, `LAG-3`, `PD-1`, TIGIT, BTLA), as.numeric)) |> 
  pivot_longer(cols = c(where(is.numeric),-Donor),
               names_to = "ICP",
               values_to = "Per_ICP") |> 
  mutate(ICP = factor(ICP, levels = c("TIM-3","BTLA","LAG-3","PD-1","CD16","CD39","TIGIT","CD161")))

#### Non Paired ###
# Df_icp_on_use_comparison <- Df_icp_on_use |>
#   filter(Purif == FALSE) |>
#   mutate(across(c(CD161, CD16, CD39, `TIM-3`, `LAG-3`, `PD-1`, TIGIT, BTLA), as.numeric)) |> 
#   pivot_longer(cols = c(where(is.numeric),-Donor),
#                names_to = "ICP",
#                values_to = "Per_ICP") |> 
#   mutate(ICP = factor(ICP, levels = c("TIM-3","BTLA","LAG-3","PD-1","CD16","CD39","TIGIT","CD161")))

#Graph
plot_icp_on_use_comparison_thaw <- ggplot(Df_icp_on_use_comparaison_paired, aes(x = Thaw, y = Per_ICP)) +
  geom_line(aes(group = Donor), alpha = 0.15) +
  facet_wrap(~ICP, strip.position = "top", nrow = 2 )+
  geom_jitter(aes(color = Thaw), width = 0, height = 0, size = 1.7, alpha = 0.75) +
  scale_color_manual(values = palette_TDG_pastel) +
  ggnewscale::new_scale_color() +
  stat_summary(aes(color = Thaw),fun = mean,geom = "crossbar", width = 0.6,
    linewidth = 1.2,
    fatten = 1) +
  scale_color_manual(values = palette_TDG) +
  stat_compare_means(method = "wilcox.test",aes(group = Thaw), hide.ns = FALSE,
                     paired = TRUE,label = "p.signif", size = 3, vjust = 1, hjust = 0.5, label.x.npc = "middle") +
  theme_blood() +
  theme(strip.text = element_text(size = 8, face = "bold"),
        legend.title = element_blank(),
        axis.title.x = element_blank())+
  labs(
    x = "Day of use",
    y = "% of CD3γ9+ cells expressing ICP \n at day of use")

#### Comparison Purif ####
Df_icp_on_use_comparaison_purif <- Df_icp_on_use |>
  mutate(across(c(CD161, CD16, CD39, `TIM-3`, `LAG-3`, `PD-1`, TIGIT, BTLA), as.numeric)) |> 
  pivot_longer(cols = c(where(is.numeric),-Donor),
               names_to = "ICP",
               values_to = "Per_ICP") |> 
  mutate(ICP = factor(ICP, levels = c("TIM-3","BTLA","LAG-3","PD-1","CD16","CD39","TIGIT","CD161"))) |> 
  mutate(Purif = if_else(Purif == TRUE, "Purified", "Non Purified"))   |> 
  mutate(Purification = if_else(Purif == "Purified", "Purified", Thaw),
         Purification = factor(Purification, levels = c("Fresh","Thaw","Purified")))   |> 
  filter(!is.na(Per_ICP)) |> 
  group_by(Donor,ICP) |>
  filter(n_distinct(Purif) == n_distinct(Df_icp_on_use$Purif)) |>
  ungroup() 

plot_icp_on_use_comparison_purification <- ggplot(Df_icp_on_use_comparaison_purif, aes(x = Purif, y = Per_ICP)) +
  geom_line(aes(group = Donor), alpha = 0.15) +
  facet_wrap(~ICP, strip.position = "top", nrow = 2 )+
  geom_jitter(aes(color = Purification), width = 0, height = 0, size = 1.7, alpha = 0.75) +
  scale_color_manual(values = palette_TDG_pastel) +
  ggnewscale::new_scale_color() +
  stat_summary(aes(color = Purif),fun = mean,geom = "crossbar", width = 0.6,
               linewidth = 1.2,
               fatten = 1) +
  scale_color_manual(values = palette_TDG) +
  stat_compare_means(method = "wilcox.test", aes(group = Purif), hide.ns = FALSE,
                     paired = TRUE, label = "p.signif", size = 3, vjust = 1, hjust = 0.5, label.x.npc = "middle") +
  theme_blood() +
  theme(strip.text = element_text(size = 8, face = "bold"),
        legend.title = element_blank(),
        axis.title.x = element_blank())+
  labs(
    x = "Day of use",
    y = "% of CD3γ9+ cells expressing ICP \n at day of use")

#### Phenotype differentiation ####
Df_pheno_differentiation <- import("./data/DF/260515_Datas_Differenciation.xlsx") |> 
  mutate(Donor = recode(Donor,"BC240112A*_1" = "BC240112A_1*",
                        "BC230331A*_1" = "BC230331A_1*",
                        "BC231020A*_1" = "BC231020A_1*",
                        "BC240322A*_1" = "BC240322A_1*",
                        "BC240209B*_1" = "BC240209B_1*")) |> 
  mutate(Thaw = if_else(str_detect(Donor, "\\*$"), "Thaw", "Fresh")) |> 
  mutate(Purification = if_else(Purified == "Yes", "Purified", Thaw),
         Purification = factor(Purification, levels = c("Fresh","Thaw","Purified"))) |> 
  pivot_longer(cols = c(`Naive CD45RA+ CCR7+`,`Central memory CD45RA- CCR7+`,
                        `Effector memory CD45RA- CCR7-`,`Terminal effector memory CD45RA+ CCR7-`),
               names_to = "Phenotype",
               values_to = "Per_phenotype") |> 
  mutate(Phenotype = factor(Phenotype, levels = c("Naive CD45RA+ CCR7+",
                                                  "Central memory CD45RA- CCR7+",
                                                  "Effector memory CD45RA- CCR7-",
                                                  "Terminal effector memory CD45RA+ CCR7-")),
         Per_phenotype = as.numeric(Per_phenotype)) |> 
  mutate(Phenotype = recode(Phenotype,"Naive CD45RA+ CCR7+" = "Naive \n CD45RA+ CCR7+",
                                           "Central memory CD45RA- CCR7+" = "Central memory \n CD45RA- CCR7+",
                                           "Effector memory CD45RA- CCR7-"= "Effector memory \n CD45RA- CCR7-",
                                           "Terminal effector memory CD45RA+ CCR7-" = "Terminal effector memory \n CD45RA+ CCR7-")) |>
  summarize(meanPurif = mean(Per_phenotype, na.rm = TRUE), .by = c(Purification, Phenotype))

plot_differentiation_tgd <- ggplot(Df_pheno_differentiation, aes(x = "", y = meanPurif, fill = Phenotype)) +
  geom_col(width = 1, color = "white", linewidth = 0.5) +
  coord_polar("y", start = 0) +
  facet_wrap(~Purification, strip.position = "top") +
  theme_blood() +
  scale_fill_manual(values = palette_diff,
                    guide = guide_legend(nrow = 2))+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 8, face = "bold"),
        legend.position = "bottom")+
  labs(fill = "Phenotype",
       title = "Differentiation of CD3γ9+ cells at day of use")
  


#### Montage ####
Figure_2 <- plot_grid(
  plot_grid(plot_Per_d1_d2_comparison,plot_icp_on_use_comparison_thaw, 
            ncol = 2, labels = c("A","B"), vjust = 1,rel_widths = c(0.55,1.45), align = "hv", axis = "bt"),
  plot_grid(plot_icp_on_use_comparison_purification,plot_differentiation_tgd,
            ncol = 2, labels = c("C","D"), rel_widths = c(1.15,0.85), axis ="bt"), # Deuxième ligne
  nrow = 2,
  axis = "bt",
  align = "hv"
)

print(Figure_2)

ggsave(
  filename = "./outputs/Figures/Figure_Lea_2.png",
  plot = Figure_2,
  device = "png",
  width = 34, # largeur A4 en cm
  height = 20, # hauteur A4 en cm
  units = "cm",
  dpi = 600,
  bg = "white"
)
