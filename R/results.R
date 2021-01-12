#load libraries
if( !require( "pacman" ) ) { install.packages( "pacman" ) }
pacman::p_load(
  tidyverse,
  here,
  ggforce
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

Du2019_end =  master_ds %>% filter(
  flask == "late"
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

Tenaillon2016_end = master_ds %>% filter(
  paper == "Tenaillon2016" &
    generation == "50000"
)

Wielgoss2016_end = master_ds %>% filter(
  paper == "Wielgoss2016" &
    generation == "10"
)

end_point_ds = rbind(
  non_mtl_wide_ds,
  Du2019_end,
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

end_point_chyper_plot = ggplot(
  data = end_point_ds, 
  aes(x = c_hyper, y = species, color = species),
  show.legend = FALSE
) +
  geom_violin() +
  geom_point() +
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position = "bottom"
  ) 
end_point_ds =end_point_ds %>% filter(species != "P_aeruginosa_PA14")
end_point_estimate_plot = ggplot(
  data = end_point_ds, 
  aes(x = species, y = estimate, color = species),
  show.legend = FALSE
) +
  geom_sina(
    stat="sina",
    position = "dodge"
  ) +
  geom_point() +
  theme(
    axis.text.x = element_text(angle = 90),
    legend.position = "bottom"
  ) 

#multiple wide - generation analysis
multiple_wide_ds = master_ds %>% 
  filter(
    func == "multiple_wide"
  ) %>% drop_na(generation)

multiple_chyper_plot = ggplot() +
  geom_point(
    data = multiple_wide_ds, 
    aes(x = generation, y = c_hyper, color = paper),
    show.legend = FALSE
  ) + theme(
    axis.text.x = element_text(angle = 90)
  ) + facet_wrap(
    .~paper,
    scales = "free"
  )+
  ylim(c(0,30))

multiple_estimate_plot = ggplot() +
  geom_point(
    data = multiple_wide_ds, 
    aes(x = generation, y = estimate, color = paper),
    show.legend = FALSE
  ) + theme(
    axis.text.x = element_text(angle = 90)
  ) + facet_wrap(
    .~paper,
    scales = "free"
  ) + ylim(c(0,1))

multiple_chyper_species = ggplot() +
  geom_point(
    data = multiple_wide_ds, 
    aes(x = generation, y = c_hyper, color = species)
  ) + theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90)
  ) + facet_wrap(
    .~paper,
    scales = "free"
  )

multiple_estimate_species = ggplot() +
  geom_point(
    data = multiple_wide_ds, 
    aes(x = generation, y = estimate, color = species)
  ) + theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90)
  ) + facet_wrap(
    .~paper,
    scales = "free"
  )

violin_multiple_chyper_species = ggplot(
  multiple_wide_ds, aes(species, c_hyper)
) + geom_violin() + geom_point()

violin_multiple_estimate_species = ggplot(
  multiple_wide_ds, aes(species, estimate)
) + geom_violin() + geom_point()

violin_multiple_chyper_paper = ggplot(
  multiple_wide_ds, aes(c_hyper, paper)
) + geom_violin() + geom_point()

violin_multiple_estimate_paper = ggplot(
  multiple_wide_ds, aes(estimate, paper)
) + geom_violin() + geom_point()