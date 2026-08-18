####Figure 1 Papier Léa TGD ####
source("./src/00_init_project.R")

##### Workflow ##### 
Workflow_png <- magick::image_read("./data/Workflow.png")

workflow_plot <- ggplot() +
  annotation_custom(
    grid::rasterGrob(as.raster(Workflow_png)),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) +
  theme_void()

##### Pheno J0####
Df_CD3g9_upon_culture <- import("./data/DF/260817_Datas_CD3g9_upon_culture.ods") |> 
  mutate(Thaw = if_else(str_detect(Donor, "\\*$"), "Thaw", "Fresh")) |> 
  pivot_longer(cols = where(is.numeric),
               names_to = "Days",
               values_to = "Per_CD3_g9") |> 
  mutate(Days = as.numeric(Days)) |> 
  filter(Thaw == "Fresh", !is.na(Per_CD3_g9) )

Df_CD3g9_upon_culture_D0 <- Df_CD3g9_upon_culture |> 
  filter(Days == "0") |> 
  mutate(Days = "Day 0")


plot_CD3g9_D0 <- ggplot(Df_CD3g9_upon_culture_D0,aes(x = Days, y = Per_CD3_g9)) +
  geom_jitter(width = 0.20,
    height = 0,
    size = 2,
    alpha = 0.75,
    color = palette_TDG_pastel["Frais"]) +
  stat_summary(fun = mean,
    geom = "crossbar",
    width = 0.6,
    linewidth = 1.2,
    fatten = 1,
    color = palette_TDG["Frais"]) +
  theme_blood() +
  theme(axis.text.x = element_blank())+
  labs(
    x = "Day 0",
    y = "% of CD3γ9+ cells in PBMC")

##### Etude Longitudinal ####
List_Donor_longitudinal <- Df_CD3g9_upon_culture |> 
  group_by(Donor) |> 
  mutate(count = n()) |> 
  filter(count > 14) |> 
  distinct(Donor) |> 
  pull(Donor)

Df_CD3g9_upon_culture_longitudinal <-  Df_CD3g9_upon_culture |> 
  filter(Donor %in% List_Donor_longitudinal ) |> 
  filter(!Days %in% c("9","12","16","23","30"))
  
plot_culture_longitudinal <- ggplot(Df_CD3g9_upon_culture_longitudinal, 
                                    aes(x = Days, y = Per_CD3_g9)) +
  geom_vline(xintercept = 17, linetype = "dashed", linewidth = 0.6)+
  geom_vline(xintercept = 24, linetype = "dashed", linewidth = 0.6)+
  geom_line(aes(group = Donor),
    linewidth = 0.7,
    alpha = 0.15,
    color = palette_TDG_pastel["Frais"]) +
  geom_point(size = 2,
    alpha = 0.6,
    color = palette_TDG_pastel["Frais"]) +
  stat_summary(aes(group = 1),
    fun = mean,
    geom = "line",
    linewidth = 1,
    color = palette_TDG["Frais"]) +
  stat_summary(aes(group = 1), fun.data = mean_se, geom = "errorbar",
    width = 0.5,
    linewidth = 0.7,
    color = palette_TDG["Frais"]) +
  theme_blood() +
  labs(x = "Days in culture",
    y = "CD3γ9+ cells (%)",
    title = "") +
  scale_x_continuous(breaks = seq(0, 30, by = 5)) +
  theme(legend.position = "none")
  

##### % d1 or d2 ####
Df_Per_d1_d2 <- import("./data/DF/260817_Datas_CD3_d1_or_d2.ods") |> 
    mutate(Donor = recode(Donor,"BC240112A*_1" = "BC240112A_1*")) |> 
    mutate(Thaw = if_else(str_detect(Donor, "\\*$"), "Thaw", "Fresh")) |> 
    pivot_longer(cols = where(is.numeric),
                 names_to = "Type_delta",
                 values_to = "Per_delta") |> 
    mutate(Type_delta = recode(Type_delta, "%CD3 d1" = "CD3+ δ1+", "%CD3 d2" = "CD3+ δ2+"))
  
Df_Per_d1_d2_fresh <-   Df_Per_d1_d2 |> 
  filter(Thaw == "Fresh", !is.na(Per_delta))
    
plot_Per_d1_d2 <- ggplot(Df_Per_d1_d2_fresh, aes(x = Type_delta, y = Per_delta)) +
    geom_jitter(width = 0.20, 
                height = 0, 
                size = 2, 
                alpha = 0.75, 
                color = palette_TDG_pastel["Frais"]) +
    stat_summary(fun = mean, 
                 geom = "crossbar", 
                 width = 0.6, 
                 linewidth = 1.2, 
                 fatten = 1, 
                 color = palette_TDG["Frais"]) +
  stat_compare_means(
    method = "wilcox.test", aes(group = Type_delta), hide.ns = FALSE, paired = TRUE,
    label = "p.signif", size = 3, vjust = 1, hjust = 0.5, label.x.npc = "middle")+
    theme_blood() +
  theme(axis.text.x = element_text(size = 8, face = "bold", color = "black")) +
    labs(x = "", y = "% of CD3γ9 cells at day of use")

### Alternative BoxPlot
# ggplot(Df_Per_d1_d2_fresh, aes(x = Type_delta, y = Per_delta)) +
#   geom_boxplot(position = position_dodge(0.8), outliers = FALSE, colour = "black", size = 0.5,
#                fill = palette_TDG_pastel["Frais"] ) +
#   geom_jitter(width = 0.20, 
#               height = 0, 
#               size = 1.8, 
#               alpha = 0.75, 
#               color = palette_TDG["Frais"]) +
#   stat_compare_means(
#     method = "wilcox.test", aes(group = Type_delta), hide.ns = FALSE, paired = TRUE,
#     label = "p.signif", size = 3, vjust = 1, hjust = 0.5, label.x.npc = "middle")+
#   theme_blood() +
#   theme(axis.text.x = element_text(size = 8, face = "bold", color = "black")) +
#   labs(x = "", y = "% of CD3γ9 cells at D17")

##### Evolution ICP #####
Df_icp_upon_culture <- import('./data/DF/260817_Datas_ICP_culture.ods') |> 
  rename(Days = "Day of culture") |> 
  pivot_longer(cols = c(where(is.numeric),-Days),
               names_to = "ICP",
               values_to = "Per_ICP") |> 
  mutate(ICP = factor(ICP, levels = c("TIM3","BTLA","LAG3","PD1","CD16","CD39","TIGIT")))

plot_icp_upon_culture <- ggplot(Df_icp_upon_culture,aes(x = Days, y = Per_ICP, color = ICP)) +
  geom_vline(xintercept = 17, linetype = "dashed", linewidth = 0.6)+
  geom_vline(xintercept = 24, linetype = "dashed", linewidth = 0.6)+
  geom_line(aes(group = ICP),linewidth = 1, alpha = 0.7) +
  geom_point(size = 2, alpha = 1) +
  theme_blood() +
  scale_color_manual(values = palette_ICP) +

  labs(x = "Days of culture",
    y = "ICP+ CD3γ9+ cells (%)",
    title = "",
    color = "") +
  scale_x_continuous(breaks = seq(0, 30, by = 5)) +
  theme(legend.position = "right")

##### Per ICP at Day of use ####



#### Montage ####
plot_CD3g9_D0 
plot_culture_longitudinal


Figure_1 <- plot_grid(
  plot_grid(workflow_plot,plot_CD3g9_D0, plot_culture_longitudinal, 
            ncol = 3, labels = c("A","B","C"), vjust = 1,rel_widths = c(1,0.5,1.5), align = "h"),
  plot_grid(plot_Per_d1_d2,plot_icp_upon_culture, 
            ncol = 2, labels = c("D","E"), rel_widths = c(0.4,1.6), align = "h"), # Deuxième ligne
  nrow = 2,
  axis = "bt",
  align = "hv"
  )

print(Figure_1)
  
  ggsave(
    filename = "./outputs/Figures/Figure_Lea_1.png",
    plot = Figure_1,
    device = "png",
    width = 25, # largeur A4 en cm
    height = 7.25, # hauteur A4 en cm
    units = "cm",
    dpi = 600,
    bg = "white"
  )
  