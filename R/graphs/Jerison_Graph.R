library(readr)
library(ggplot2)
OT <-read_csv("Jersion2017_Analysis_OT.csv")
HT<-read_csv("Jersion2017_Analysis_HT.csv")

OT <-data.frame(OT$Founder, OT$c_hyper) 
HT <- data.frame(HT$Founder, HT$c_hyper)
OT.1 <- cbind(OT, rep("OT",32))
HT.1 <- cbind(HT, rep("HT", 34))
colnames(OT.1) <- c("Founder","C_hyper","Temperature")
colnames(HT.1) <- c("Founder","C_hyper","Temperature")
Data <- rbind(OT.1, HT.1)

graph_Jerison <- 
  ggplot(Data, aes(x=Founder , y = C_hyper , colour= Temperature)) +
  geom_point(size = 3, alpha=0.5) + 
  theme(axis.text=element_text(size=13, angle = 90, hjust = 1), axis.title=element_text(size=16,face="bold")) + 
  xlab("Founder") + 
  ylab("C_hyper value")+coord_cartesian(ylim = c(-5, 81)) 


png("/Users/Margot/Desktop/Co-op/data-out/graph_Jerison.png", height=10, width = 10, units="in", res=300)
graph_Jerison
dev.off()
system("open /Users/Margot/Desktop/Co-op/data-out/graph_Jerison.png")
