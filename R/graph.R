library(readr)
Lang <-read_csv("data-out/Lang2014_Analysis.csv")
Sherlock2013 <- read_csv("data-out/Sherlock2013_Analysis.csv")
Sherlock2019 <- read_csv("data-out/Sherlock2019_Analysis.csv")
Kacar <- read_csv("Kacar2017_Analysis.csv")
X.label <- c(0,7, 133,140,196,240,266,322,335,385,415,448,505,585,665,745,825,910,1000)
generations <- c("7", "70", "133", "196", "266", "322", "385", "448")
Sherlock2013$generations<- generations

Lang <-data.frame(Lang$`as.numeric(generations)`, Lang$c_hyper) 
Sherlock2013 <- data.frame(Sherlock2013$generations, Sherlock2013$c_hyper)
#Sherlock2019 <- data.frame(Sherlock2019$generations, Sherlock2019$c_hyper)
#Kacar <-data.frame(Kacar$generations, Kacar$c_hyper)
#Sherlock2019.1 <- cbind(Sherlock2019, rep("Sherlock2019", 11))
Lang.1 <- cbind(Lang, rep("Lang",12))
Sherlock2013.1 <- cbind(Sherlock2013, rep("Sherlock2013", 8))
#Kacar.1 <- cbind(Kacar, rep("Kacar", 4))
#colnames(Kacar.1) <- c("Generation","C_hyper","Paper")
colnames(Sherlock2013.1) <- c("Generation","C_hyper","Paper")
#colnames(Sherlock2019.1) <- c("Generation","C_hyper","Paper")
colnames(Lang.1) <- c("Generation","C_hyper","Paper")
Data <- rbind(Lang.1, Sherlock2013.1)

graph2 <- ggplot(Data, aes(x=Generation , y = C_hyper , colour= Paper)) +geom_point(size = 3) + theme(axis.text=element_text(size=13), axis.title=element_text(size=16,face="bold")) + xlab("Time(generations)") + ylab("C_hyper value")+ expand_limits(y=c(0, 30))

setwd("/Users/Margot/Desktop/Co-op/data-out/")
png("graph2.png",height=10, width = 10, units="in", res=300)
graph2
dev.off()

  