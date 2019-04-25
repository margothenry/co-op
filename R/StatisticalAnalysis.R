library(here)
library(tidyverse)
MetaData<-read.csv(here("data-out","data-out_info", "MetaData.csv"))

#paired t-test
t.test(x = MetaData$c.score, y = MetaData$generation.test, paired = TRUE)
t.test(y = MetaData$c.score, x = MetaData$generation.test, paired = TRUE)
# Paired t-test
# 
# data:  MetaData$generation.test and MetaData$c.score
# t = 3.043, df = 138, p-value = 0.002805
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   556.0821 2619.6089
# sample estimates:
#   mean of the differences 
# 1587.845 
MetaData %>% 
  ggplot(aes(y= c.score, x= generation.test)) + geom_violin(fill = "lightyellow") + geom_point()

#Two samepl T-test
t.test(y = MetaData$c.score, x = MetaData$generation.test, var.equal = TRUE)
# Two Sample t-test
# 
# data:  MetaData$generation.test and MetaData$c.score
# t = 3.0647, df = 278, p-value = 0.002393
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   567.9036 2607.6011
# sample estimates:
#   mean of x mean of y 
# 1598.3309   10.5786 

#Welch's t-test
t.test(y = MetaData$c.score, x = MetaData$generation.test, var.equal = FALSE)
# Welch Two Sample t-test
# 
# data:  MetaData$generation.test and MetaData$c.score
# t = 3.0428, df = 138, p-value = 0.002807
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   555.9903 2619.5143
# sample estimates:
#   mean of x mean of y 
# 1598.3309   10.5786 

#One sample t-test
t.test(MetaData$generation.test)
# One Sample t-test
# 
# data:  MetaData$generation.test
# t = 3.0631, df = 138, p-value = 0.002635
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#   566.5714 2630.0905
# sample estimates:
#   mean of x 
# 1598.331 

t.test(MetaData$c.score)
# One Sample t-test
# 
# data:  MetaData$c.score
# t = 9.1794, df = 140, p-value = 5.152e-16
# alternative hypothesis: true mean is not equal to 0
# 95 percent confidence interval:
#   8.300179 12.857023
# sample estimates:
#   mean of x 
# 10.5786 

#wilcoxon matched pairs
wilcox.test(y = MetaData$c.score, x = MetaData$generation.test, paired = TRUE)
# Wilcoxon signed rank test with continuity correction
# 
# data:  MetaData$generation.test and MetaData$c.score
# V = 9557, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

wilcox.test(y = MetaData$c.score, x = MetaData$generation.test)
# Wilcoxon rank sum test with continuity correction
# 
# data:  MetaData$generation.test and MetaData$c.score
# W = 18668, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

# correlation
cor.test(y = MetaData$c.score, x = MetaData$generation.test, method = "pearson")
# Pearson's product-moment correlation
# 
# data:  MetaData$generation.test and MetaData$c.score
# t = -0.006296, df = 137, p-value = 0.995
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.1670238  0.1659778
# sample estimates:
#           cor 
# -0.0005379068 

cor.test(y = MetaData$c.score, x = MetaData$generation.test, method = "spearman")
# Spearman's rank correlation rho
# 
# data:  MetaData$generation.test and MetaData$c.score
# S = 452070, p-value = 0.9067
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#         rho 
# -0.01003407 

cor.test(y = MetaData$c.score, x = MetaData$generation.test, method = "kendall")
# Kendall's rank correlation tau
# 
# data:  MetaData$generation.test and MetaData$c.score
# z = -0.14257, p-value = 0.8866
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#          tau 
# -0.009118304 

#goodness of fit
chisq.test(MetaData$c.score)
# Chi-squared test for given probabilities
# 
# data:  MetaData$c.score
# X-squared = 2478.3, df = 140, p-value < 2.2e-16
