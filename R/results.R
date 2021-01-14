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
    generation == "50000"
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

end_point_chyper_plot = ggplot(
  data = end_point_ds, 
  aes(x = species, y = c_hyper, color = species),
  show.legend = FALSE
) +
  geom_sina(
    stat = "sina",
    position = "dodge"
  ) +
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
    aes(x = generation, y = c_hyper, color = species),
    show.legend = TRUE
  ) + theme(
    axis.text.x = element_text(angle = 90),
    legend.position = "bottom"
  ) + facet_wrap(
    .~paper,
    scales = "free"
  )+
  ylim(c(0,36))


violin_multiple_chyper_species = ggplot(
  multiple_wide_ds, aes(species, c_hyper)
) + geom_sina(
  stat = "sina",
  position = "dodge"
) 

violin_multiple_estimate_species = ggplot(
  multiple_wide_ds, aes(species, estimate)
) +  geom_sina(
  stat = "sina",
  position = "dodge"
) 

# generation by chyper by species
