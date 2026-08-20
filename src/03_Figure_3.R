####Figure 3 Papier Léa TGD ####
source("./src/00_init_project.R")


##### Workflow ##### 
Workflow_malc_png <- magick::image_read("./data/Workflow_malc_tgd.png")

workflow_malc_plot <- ggplot() +
  annotation_custom(
    grid::rasterGrob(as.raster(Workflow_malc_png)),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) +
  theme_void()


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

# plot_depletion_UT_thaw_boxplot <- ggplot(Df_depletion_UT_thaw, 
#                                  aes(x = Purification, y = UTplusT_normalized_depletion_by_UTminusT)) +
#   geom_boxplot(aes(fill = Purification),width = 0.6, outliers = FALSE,
#                linewidth = 0.7,
#                color = "black") +
#   scale_fill_manual(values = palette_TDG) +
#   ggnewscale::new_scale_color() +
# geom_beeswarm(aes(color = Purification), cex = 3, size = 1.7, alpha = 0.75)+
#   scale_color_manual(values = palette_TDG_pastel) +
#   stat_compare_means(method = "wilcox.test",aes(group = Purification), hide.ns = FALSE,
#                      paired = FALSE,label = "p.signif", size = 3, vjust = 1, hjust = 0.5, label.x.npc = "middle") +
#   theme_blood() +
#   theme(strip.text = element_text(size = 8, face = "bold"),
#         legend.position =  "none") +
#   labs(x = "", y = "Normalized depletion \n by UT- T cells")

plot_depletion_UT_thaw <- ggplot(Df_depletion_UT_thaw, 
                                 aes(x = Purification, y = UTplusT_normalized_depletion_by_UTminusT)) +
  stat_summary(aes(color = Thaw),fun = mean,geom = "crossbar", width = 0.6,
               linewidth = 1.2,
               fatten = 1) +
  scale_color_manual(values = palette_TDG) +
  ggnewscale::new_scale_color() +
  geom_beeswarm(aes(color = Purification), cex = 4, size = 1.7, alpha = 0.75)+
  scale_color_manual(values = palette_TDG_pastel) +
  stat_compare_means(method = "wilcox.test",aes(group = Purification), hide.ns = FALSE,
                     paired = FALSE,label = "p.signif", size = 3, vjust = 1, hjust = 0.5, label.x.npc = "middle") +
  theme_blood() +
  theme(strip.text = element_text(size = 8, face = "bold"),
        legend.position =  "none") +
  labs(x = "", y = "Normalized depletion \n by UT- T cells")

##### Depletion Purif #####
plot_depletion_UT_Purif <- ggplot(Df_depletion_UT, aes(x = Purif, y = UTplusT_normalized_depletion_by_UTminusT)) +
  stat_summary(aes(color = Purif),fun = mean,geom = "crossbar", width = 0.6,
               linewidth = 1.2,
               fatten = 1) +
  scale_color_manual(values = palette_TDG) +
  ggnewscale::new_scale_color() +
  geom_beeswarm(aes(color = Purification), cex = 4, size = 1.7, alpha = 0.75)+
  scale_color_manual(values = palette_TDG_pastel) +
  stat_compare_means(method = "wilcox.test",aes(group = Purif), hide.ns = FALSE,
                     paired = FALSE,label = "p.signif", size = 3, vjust = 1, hjust = 0.5, label.x.npc = "middle") +
  theme_blood() +
  theme(strip.text = element_text(size = 8, face = "bold"),
        legend.position =  "none") +
  labs(x = "", y = "Normalized depletion \n by UT- T cells")

# plot_depletion_UT_Purif_boxplot <- ggplot(Df_depletion_UT, aes(x = Purif, y = UTplusT_normalized_depletion_by_UTminusT)) +
#   geom_boxplot(aes(fill = Purif),width = 0.6, outliers = FALSE,
#                linewidth = 0.7,
#                color = "black") +
#   scale_fill_manual(values = palette_TDG) +
#   ggnewscale::new_scale_color() +
# geom_beeswarm(aes(color = Purification), cex = 4, size = 1.7, alpha = 0.75)+
#   scale_color_manual(values = palette_TDG_pastel) +
#   stat_compare_means(method = "wilcox.test",aes(group = Purif), hide.ns = FALSE,
#                      paired = FALSE,label = "p.signif", size = 3, vjust = 1, hjust = 0.5, label.x.npc = "middle") +
#   theme_blood() +
#   theme(strip.text = element_text(size = 8, face = "bold"),
#         legend.position =  "none") +
#   labs(x = "", y = "Normalized depletion \n by UT- T cells")

#### GrzmB / CD107a ####
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
  rename( "Tγδ only"= "UT_T_only",
          "Tγδ only + GA101"= "GA_Tonly",
          "Tγδ + RL" = "UTplusT",
          "Tγδ + RL + GA101" = "GAplusT")

p <- ggplot(Df_cd107a, aes(x = Purification, y = `Tγδ + RL`)) +
  stat_summary(aes(color = Purif),fun = mean,geom = "crossbar", width = 0.6,
               linewidth = 1.2,
               fatten = 1) +
  scale_color_manual(values = palette_TDG) +
  ggnewscale::new_scale_color() +
  geom_beeswarm(aes(color = Purification), cex = 4, size = 1.7, alpha = 0.75)+
  scale_color_manual(values = palette_TDG_pastel) +
  stat_compare_means(method = "wilcox.test",aes(group = Purification), hide.ns = FALSE,
                     paired = FALSE,label = "p.signif", size = 3, vjust = 1, hjust = 0.5, label.x.npc = "middle")+
  labs(color = "",
       fill ="")+
  theme_blood()

legend_fig <- get_legend(p)
plot_grid(legend_fig)
##### Non Purified ######
Df_cd107a_non_purified <- Df_cd107a |> 
  filter(Purified == "No", !is.na(`Tγδ only`)) |> 
  pivot_longer(cols = c(where(is.numeric)),
               names_to = "Condition",
               values_to = "Per_CD107a") |> 
  filter(Condition %in% c("Tγδ only", "Tγδ + RL")) |> 
  mutate(Condition = factor(Condition, levels = c("Tγδ only", "Tγδ + RL")))

# plot_CD107a_non_purified_boxplot <- ggplot(Df_cd107a_non_purified, aes(x = Condition, y = Per_CD107a)) +
#   geom_line(aes(group = Donor), alpha = 0.15) +
#   geom_boxplot(aes(fill = Purif),width = 0.6, outliers = FALSE,
#                linewidth = 0.7,
#                color = "black") +
#   scale_fill_manual(values = palette_TDG) +
#   ggnewscale::new_scale_color() +
#   geom_beeswarm(aes(color = Purification), cex = 4, size = 1.7, alpha = 0.75)+
#   scale_color_manual(values = palette_TDG_pastel) +
#   stat_compare_means(method = "wilcox.test",aes(group = Condition), hide.ns = FALSE,
#                      paired = TRUE,label = "p.signif", size = 3, vjust = 1, hjust = 0.5, label.x.npc = "middle") +
#   theme_blood() +
#   theme(strip.text = element_text(size = 8, face = "bold"),
#         legend.position =  "none") +
#   labs(x = "", y = "% CD107a")

plot_CD107a_non_purified <- ggplot(Df_cd107a_non_purified, aes(x = Condition, y = Per_CD107a)) +
  geom_line(aes(group = Donor), alpha = 0.15) +
  stat_summary(aes(color = Purif),fun = mean,geom = "crossbar", width = 0.6,
               linewidth = 1.2,
               fatten = 1) +
  scale_color_manual(values = palette_TDG) +
  ggnewscale::new_scale_color() +
  geom_beeswarm(aes(color = Purification), cex = 2.5, size = 1.7, alpha = 0.75)+
  scale_color_manual(values = palette_TDG_pastel) +
  stat_compare_means(method = "wilcox.test",aes(group = Condition), hide.ns = FALSE,
                     paired = TRUE,label = "p.signif", size = 3, vjust = 1, hjust = 0.5, label.x.npc = "middle") +
  theme_blood() +
  ylim(0,50)+
  theme(strip.text = element_text(size = 8, face = "bold"),
        legend.position =  "none") +
  labs(x = "", y = "% CD107a")

##### Purified ######
Df_cd107a_purified <- Df_cd107a |> 
  filter(Purified == "Yes", !is.na(`Tγδ only`)) |> 
  pivot_longer(cols = c(where(is.numeric)),
               names_to = "Condition",
               values_to = "Per_CD107a") |> 
  filter(Condition %in% c("Tγδ only", "Tγδ + RL")) |> 
  mutate(Condition = factor(Condition, levels = c("Tγδ only", "Tγδ + RL")))

# plot_CD107a_purified_boxplot <- ggplot(Df_cd107a_purified, aes(x = Condition, y = Per_CD107a)) +
#   geom_line(aes(group = Donor), alpha = 0.15) +
#   geom_boxplot(aes(fill = Purif),width = 0.6, outliers = FALSE,
#                linewidth = 0.7,
#                color = "black") +
#   scale_fill_manual(values = palette_TDG) +
#   ggnewscale::new_scale_color() +
#   geom_beeswarm(aes(color = Purification), cex = 1, size = 1.7, alpha = 0.75)+
#   scale_color_manual(values = palette_TDG_pastel) +
#   stat_compare_means(method = "wilcox.test",aes(group = Condition), hide.ns = FALSE,
#                      paired = TRUE,label = "p.signif", size = 3, vjust = 1, hjust = 0.5, label.x.npc = "middle") +
#   theme_blood() +
#   theme(strip.text = element_text(size = 8, face = "bold"),
#         legend.position =  "none") +
#   labs(x = "", y = "% CD107a")

plot_CD107a_purified <- ggplot(Df_cd107a_purified, aes(x = Condition, y = Per_CD107a)) +
  geom_line(aes(group = Donor), alpha = 0.15) +
  stat_summary(aes(color = Purif),fun = mean,geom = "crossbar", width = 0.6,
               linewidth = 1.2,
               fatten = 1) +
  scale_color_manual(values = palette_TDG) +
  ggnewscale::new_scale_color() +
  geom_beeswarm(aes(color = Purification), cex = 2.5, size = 1.7, alpha = 0.75)+
  scale_color_manual(values = palette_TDG_pastel) +
  stat_compare_means(method = "wilcox.test",aes(group = Condition), hide.ns = FALSE,
                     paired = TRUE,label = "p.signif", size = 3, vjust = 1, hjust = 0.5, label.x.npc = "middle") +
  theme_blood() +
  ylim(0,50)+
  theme(strip.text = element_text(size = 8, face = "bold"),
        legend.position =  "none") +
  labs(x = "", y = "% CD107a on CD3+ Tγδ")

#### Imagerie ####
Figure_imagerie_png <- magick::image_read("./data/Figure_imagerie.png")

Figure_imagerie_plot <- ggplot() +
  annotation_custom(
    grid::rasterGrob(as.raster(Figure_imagerie_png)),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) +
  theme_void()

#### Montage ####
Figure_3 <- plot_grid(
  plot_grid(workflow_malc_plot,Figure_imagerie_plot,
            ncol = 2, labels = c("A"), vjust = 1.1,rel_widths = c(0.95,1.05), align = "hv", axis = "bt"),
  plot_grid(plot_depletion_UT_thaw, plot_depletion_UT_Purif,plot_CD107a_non_purified,plot_CD107a_purified,legend_fig,
            ncol = 5, labels = c("E","F","G","H",""), rel_widths = c(0.9,0.9,0.9,0.9,0.4),
            vjust = 0.8), # Deuxième ligne
  nrow = 2,
  rel_heights = c(1,0.5),
  axis = "bt",
  align = "hv"
)

print(Figure_3)

ggsave(
  filename = "./outputs/Figures/Figure_Lea_3.png",
  plot = Figure_3,
  device = "png",
  width = 34, # largeur A4 en cm
  height = 20, # hauteur A4 en cm
  units = "cm",
  dpi = 600,
  bg = "white"
)



