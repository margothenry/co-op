library(tidyverse)
#Testing if the median in E.coli is different then the median in Sac for their c_yper value.
#SOOOO our data can not use anova since it doesnt meet the assumptions of normality or equal variances. so we can use the wilcox test. 
data <- read_csv("~/Desktop/stat 4690/STAT4690/end_point_generation_analysis.csv")
 
# data <- data %>% select(species, c_hyper)
# data <- data %>% group_by(species) %>% mutate(X= row_number()) %>% 
#   spread(species, c_hyper) %>% select(-X) 

summary = data %>% group_by( species) %>%
  summarise(
    count = n(),
    median = median(c_hyper, na.rm = TRUE),
    IQR = IQR(c_hyper, na.rm = TRUE)
  )

library("ggpubr")
pic <- ggboxplot(data, x = "species", y = "c_hyper", 
          color = "species", palette = c("#00AFBB", "#E7B800"),
          ylab = "C_hyper", xlab = "Species")

test <- wilcox.test(c_hyper ~ species, data = data,
                   exact = FALSE)
test
#reject Ho big time
#We can conclude that E.coli's median c_hyper is significantly different from Sac's median c_hyper with a p-value = 3.957e-07, which correcsponds with the graph above.