#load libraries
if( !require( "pacman" ) ) { install.packages( "pacman" ) }
pacman::p_load(
  tidyverse,
  here,
  ggforce,
  ggtext
)

master_ds <- read.csv(
  here("data_out", "master_analyses.csv"),
  encoding = "UTF-8"
)

#end point analysis
non_mtl_wide_ds = master_ds %>% 
  filter(
    func != "multiple_wide"
  )

Flynn2014_end = master_ds %>% filter(
  paper == "Flynn2014" &
    generation == "540"
)

Kacar2017_end = master_ds %>% filter(
  paper == "Kacar2017" &
    generation == "2000"
)

Keane2014_end = master_ds %>% filter(
  paper == "Keane2014" &
    generation == "2200"
)

Lang2013_end = master_ds %>% filter(
  paper == "Lang2013" &
    generation == "1000"
)

Morgenthaler2019_end = master_ds %>% filter(
  paper == "Morgenthaler2019" &
    day == "50"
)

Sandberg2016_end = master_ds %>% filter(
  paper == "Sandberg2016" &
    flask == "133"
)

Sherlock2013_end = master_ds %>% filter(
  paper == "Sherlock2013" &
    generation == "448"
)

Sherlock2019_end = master_ds %>% filter(
  paper == "Sherlock2019" &
    generation == "500"
)

Tenaillon2016_end = master_ds %>% filter(
  paper == "Tenaillon2016" &
    generation == "2000"
)

Wielgoss2016_end = master_ds %>% filter(
  paper == "Wielgoss2016" &
    day == "10"
)

Avrani2017_end = master_ds %>% filter(
  paper == "Avrani2017" &
    day == "127"
)

end_point_ds = rbind(
  non_mtl_wide_ds,
  Avrani2017_end,
  Sherlock2019_end,
  Flynn2014_end,
  Kacar2017_end,
  Keane2014_end,
  Lang2013_end,
  Morgenthaler2019_end,
  Sandberg2016_end,
  Sherlock2013_end,
  Tenaillon2016_end,
  Wielgoss2016_end
)


end_point_ds$kingdom[end_point_ds$species == "Sac"] <- "fungi"
end_point_ds$kingdom[end_point_ds$species != "Sac"] <- "bacteria"

lm_fit <- lm(c_hyper ~ generation + kingdom , data = end_point_ds)

end_point_chyper_plot = ggplot(
  data = end_point_ds, 
  aes(x = kingdom, y = c_hyper),
  show.legend = FALSE
) +
  geom_boxplot(
    outlier.shape = NA,
    varwidth = 0.2
    ) +
  geom_jitter(
    aes(
      color = species),
    alpha = 0.5,
    show.legend = FALSE,
    position = position_jitter(width = 0.2, seed = 0),
    size = 2
    ) +
  scale_color_manual(values = c("cyan4","darkorange","purple")) +
ylab("Genetic repeatability") +
  theme_bw()

gen_end_point_ds = end_point_ds
#############
#italic and change names and send
#plot again with teniallon gen ~ 2000
#time series with linear regression line
gen_end_point_ds$species = factor(
  gen_end_point_ds$species,
  levels = c(
    "Ecoli_K12",
    "Sac",
    "P_aeruginosa_PA14"
  ),
  labels = c(
    "E. coli",
    "S. cerevisiae",
    "P. aeruginosa"
  )

) 
gen_end_point_chyper_plot = ggplot() +
  geom_point(
    data = gen_end_point_ds, 
    aes(x = generation, y = c_hyper, color = species),
    size = 2,
    alpha = 0.7
  ) +
scale_color_manual(
  values = c("cyan4","darkorange","purple"),
  labels = c(
    expression(italic("E. coli")),
    expression(italic("S. cerevisiae")),
    expression(italic("P. aeruginosa")))
  ) +
  ylab("Genetic repeatability") +
  xlab("Generation") +
  theme_bw() +
  theme(
    legend.text.align = 0
  ) 


gen_end_point_chyper_plot = ggplot() +
  geom_point(
    data = gen_end_point_ds, 
    aes(x = generation, y = c_hyper, color = species),
    size = 2,
    alpha = 0.7
  ) +
  scale_color_manual(
    values = c("cyan4","darkorange","purple"),
    labels = c(
      expression(italic("E. coli")),
      expression(italic("S. cerevisiae")),
      expression(italic("P. aeruginosa")))
  ) +
  ylab("Genetic repeatability") +
  xlab("Generation") +
  theme_bw() +
  theme(
    legend.text.align = 0
  ) 


#multiple wide - generation analysis
multiple_wide_ds = master_ds %>% 
  filter(
    func == "multiple_wide"
  ) %>% drop_na(generation)


ten_chyper_plot = ggplot() +
  geom_point(
    data = multiple_wide_ds, 
    aes(x = generation, y = c_hyper, colour = species),
    show.legend = FALSE,
    size = 2
  ) + theme(
    axis.text.x = element_text(angle = 90),
    legend.position = "bottom"
  ) + 
  ylim(c(0,30)) +
  scale_color_manual(values = c("cyan4")) +
  ylab("Genetic repeatability") +
  theme_bw() + 
  geom_abline(slope = coef(lm_fit)[[2]], intercept = coef(lm_fit)[[1]])+
  geom_smooth(
    method = "lm"
  )

gen_chyper_plot = ggplot() +
  geom_point(
    data = multiple_wide_ds, 
    aes(x = generation, y = c_hyper, colour = species),
    show.legend = FALSE,
    size = 2
  ) + theme(
    axis.text.x = element_text(angle = 90),
    legend.position = "bottom"
  ) + 
  ylim(c(0,37)) +
  scale_color_manual(values = c("darkorange","purple","cyan4")) +
  ylab("Genetic repeatability") +
  theme_bw()

#lm generations
# lm_fit = lm(c_hyper ~ generation , data = multiple_wide_ds)

multiple_chyper_plot = ggplot() +
  geom_point(
    data = multiple_wide_ds, 
    aes(x = generation, y = c_hyper, color = species)
  ) + stat_smooth(method="lm") +
  facet_wrap(
    .~paper,
    scales = "free",
    nrow = 2,
    ncol = 4
  )+
  # geom_abline(slope = coef(lm_fit)[[2]], intercept = coef(lm_fit)[[1]])+
  ylim(c(0,37)) +
  scale_color_manual(
    values = c("cyan4","darkorange","purple"),
    labels = c(
      expression(italic("E. coli")),
      expression(italic("S. cerevisiae")),
      expression(italic("P. aeruginosa")))
  ) +
  ylab("Genetic repeatability") +
  xlab("Generation") +
  theme_bw() +
  theme(
    legend.text.align = 0,
    axis.text.x = element_text(angle = 90),
    legend.position = "bottom"
  ) 
