library(here)
library(dgconstraint)
library(tidyverse)
library(R.utils)
sourceDirectory(here("R", "functions"))

############################
#multiple generations
############################

### Tenaillon2016
# This paper originally had two clones sequenced for each timepoint. Collapse the two columns and save the resulting dataset.
Tenaillon2016 <- read_csv(here("data-in", "Tenaillon2016.csv"))
Tenaillon2016 <- Tenaillon2016 %>% 
  transmute(Gene = Tenaillon2016$Gene, Population = Tenaillon2016$Population, Details = Tenaillon2016$Details, `500` = `500 I1 R1`+`500 I2 R1`, `1000` = `1000 I1 R1`+`1000 I2 R1`, `1500` =`1500 I1 R1`+`1500 I2 R1`, `2000` = `2000 I1 R1`+`2000 I2 R1`, `5000` = `5000 I1 R1`+`5000 I2 R1`, `10000` = `10000 I1 R1`+`10000 I2 R1`, `15000` = `15000 I1 R1`+`15000 I2 R1`, `20000` = `20000 I1 R1`+`20000 I2 R1`, `30000` = `30000 I1 R1`+`30000 I2 R1`, `40000` = `40000 I1 R1`+`40000 I2 R1`, `50000` = `50000 I1 R1`+`50000 I2 R1`) %>% replace(is.na(.), 0) 
write_csv(Tenaillon2016, "data-in/Tenaillon2016.csv")

multigen_c_hyper("Tenaillon2016-original", "David minimal medium", "Ecoli_K12", c("500", '1000', '1500', '2000', '5000', '10000', '20000', '30000', '40000','50000'))

### Lang2014
#filter out intergenics before running the code
Lang2014 <-read_csv(here("data-in", "Lang2014.csv"))
Lang2014 <- Lang2014 %>%
  filter(Gene != "Intergenic")

multigen_c_hyper("Lang2014","YPD", "Sac",c("0","140","240","335",'415','505','585','665','745','825','910','1000'))

### Sherlock2013
multigen_c_hyper("Sherlock2013", "YPD","Sac" ,c("7","70", "133","196","266", "322","385","448"))

### Sherlock2013
multigen_c_hyper("Sherlock2019","Davids minimal medium" ,"Ecoli_K12",c("0","50","100","150","200" ,"250","300","350","400","450","500"))

### Kacar2017
#filter out intergenics
Kacar2017 <-read_csv(here("data-in", "Kacar2017.csv"))
Kacar2017 <- Kacar2017 %>% 
            filter(`Nature of the change` != "Intergenic")

multigen_c_hyper("Kacar2017","Minimal glucose medium","Ecoli_K12", c("500", "1000", "1500", "2000"))

###Tonoyan_2019
#filter out intergenics
Tonoyan2019 <-read_csv(here("data-in", "Tonoyan2019.csv"))
Tonoyan2019 <-Tonoyan2019  %>% 
              filter(Gene != "Intergenic")

multigen_c_hyper("Tonoyan_2019", "LB medium","Ecoli_ATCC",c("1", "14", "20"))

### Wielgoss2016a
#filter out intergenics
Wielgoss2016-mucoid <-read_csv(here("data-in", "Wielgoss2016-mucoid.csv"))
Wielgoss2016-mucoid<- Wielgoss2016-mucoid%>% 
                filter(Gene != "Intergenic")

multigen_c_hyper("Wielgoss2016-mucoid", "lysogeny broth", "Ecoli_K12",c("0", "10"))

### Wielgoss2016-nonMucois
#filter out intergenics
Wielgoss2016-nonMucoid <-read_csv(here("data-in", "Wielgoss2016-nonMucoid.csv"))
Wielgoss2016-nonMucoid<- Wielgoss2016-nonMucoid%>% 
  filter(Gene != "Intergenic")

multigen_c_hyper("Wielgoss2016-nonMucoid","lysogeny broth", "Ecoli_K12", c("0","2", "10"))

### Sandberg2016
multigen_c_hyper("Sandberg2016", "Davids minimal medium", "Ecoli_K12", c("Flask 23", "Flask 58", "Flask 133"))

#####################################
# multiple selective pressures
####################################

### Hong2014
## filter out unwanted mutations
Hong2014 <-read_csv(here("data-in", "Hong2014.csv"))
Hong2014 <-Hong2014 %>% 
  filter( type != 'DIP', type !='WCA', type != 'IND', type != 'AMP')

multipressure_c_hyper("Hong2014",c("Ammonium", "Arginine", "Urea", "Allantoin"), "Sac", c("Ammonium", "Arginine", "Urea", "Allantoin"))

### Jersion2017
#for this data set there were multiple different ways to run an analysis and so each of the following code chunks will results in a different anlysis.
#interested only in the populations or clones in "OT"
Jerison2017 <- read_csv(here("data-in", "Jerison2017.csv"))
Jerison2017<-Jerison2017 %>% 
  filter(Gene != "Intergenic") %>% 
  filter(selective_pressure == "OT")
Founder <- unique(data$Founder)

#interested only in the populations or clones in "OT"
Jerison2017 <- read_csv(here("data-in", "Jerison2017.csv"))
Jerison2017<-Jerison2017 %>% 
  filter(Gene != "Intergenic") %>% 
  filter(selective_pressure == "HT")
Founder <- unique(data$Founder)

#interested in the difference in c-hyper between the two selectve pressures.
Jerison2017 <- read_csv(here("data-in", "Jerison2017.csv"))
Jerison2017<-Jerison2017 %>% 
  filter(Gene != "Intergenic")
Founder <- unique(data$Founder)

multipressure_c_hyper("Jersion2017",c("OT", "HT"),"Sac",c("OT", "HT"))

### Payen2016
#There are a few ways to run this analysis
#filter out the clones and analyze between the selective pressure whether the ploidy
Payen2016 <- read_csv(here("data-in", "Payen2016.csv"))
Payen2016 <- Payen2016 %>% 
            filter(frequency != "clone")%>% 
            drop_na(frequency) 

#analyze between the selective pressures for only the diploid genes
Payen2016 <- read_csv(here("data-in", "Payen2016.csv"))
Payen2016 <- Payen2016 %>% 
              filter(Ploidy == "diploid") %>% 
            filter(frequency != "clone")%>% 
            drop_na(frequency) 

#analyze between the selective pressures for only the haploid genes 
Payen2016 <- read_csv(here("data-in", "Payen2016.csv"))
Payen2016 <- Payen2016 %>% 
            filter(Ploidy == "diploid") %>% 
            filter(frequency != "clone")%>% 
            drop_na(frequency) 

multipressure_c_hyper("Payen2016",c("phosphate", "sulfate", "glucose"), "Sac", c("phosphate", "sulfate", "glucose"))
##############################
# Single Generation Matrix
#############################

### Wannier2018
# This paper originally had two clones sequenced for each population. Collapse the two columns and save the resulting dataset. And filter out the intergenic
Wannier2018 <- read_csv(here("data-in", "Wannier2018.csv"))
Wannier2018<- Wannier2018 %>% transmute(Gene = Wannier2018$Gene, Details= Wannier2018$Details, A1= `A1 F1 I1 R1`+ `A1 F1 I2 R1`, A2 = `A2 F1 I1 R1`+ `A2 F1 I2 R1`, A3 = `A3 F1 I1 R1`+ `A3 F1 I2 R1`, A4 = `A4 F1 I1 R1`+ `A4 F1 I2 R1`, A5 = `A5 F1 I1 R1`+ `A5 F1 I2 R1`, A6 = `A6 F1 I1 R1`+ `A6 F1 I2 R1`, A7 = `A7 F1 I1 R1`+ `A7 F1 I2 R1`, A8 = `A8 F1 I1 R1`+ `A8 F1 I2 R1`, A9 = `A9 F1 I1 R1`+ `A9 F1 I2 R1`, A10 = `A10 F1 I1 R1`+ `A10 F1 I2 R1`,A11 = `A1 F1 I1 R1`+ `A11 F1 I2 R1`, A12 = `A12 F1 I1 R1`+ `A12 F1 I2 R1`, A13 = `A13 F1 I1 R1`+ `A13 F1 I2 R1`, A14 = `A14 F1 I1 R1`+ `A14 F1 I2 R1`) %>% replace(is.na(.), 0) 
population <- c("A1" , "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12", "A13", "A14")
Wannier2018 <-Wannier2018 %>% filter(Details != "intergenic") 
Wannier2018.1<-Wannier2018 %>% select(population)
Wannier2018.matrix<- as.matrix(Wannier2018.1)
Wannier2018.matrix[Wannier2018.matrix > 0] <-1
Wannier201<- cbind(as.data.frame(Wannier2018.matrix), Wannier2018$Gene)

            
singlematrix_c_hyper("Wannier2018","Glucose Mininal Medium", "Ecoli_K12", c("A1" , "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11", "A12", "A13", "A14"))

### Long2017
singlematrix_c_hyper("Long2017", "Glucose minimal media","Ecoli_K12", c("ALE1", "ALE2", "ALE3", "ALE4", "ALE5", "ALE6", "ALE7", "ALE8", "ALE9", "ALE10"))

### Sandberg2014
singlematrix_c_hyper("Sandberg2014","Glucose minimal media", "Ecoli_K12",c("ALE1", "ALE2", "ALE3", "ALE4", "ALE5", "ALE6", "ALE7", "ALE8", "ALE9", "ALE10"))

### Creamer2016
#filter out the intergenics
Creamer2016 <-read_csv(here("data-in", "Creamer2016.csv"))
Creamer2016 <-Creamer2016 %>% 
              filter(Annotation!= "intergenic")

singlematrix_c_hyper( "Creamer2016","LBK medium","Ecoli_K12",c("K0001","K0011","K0002","K0003","K0014","K0015","KB026","K0022","K0023","K0006","K0007","K0010", "K0020","K0019","K0031","K0030"))

################################
#Single Generation NonMatrix
################################

### Deatherage2017
#This paper originally had three clones sequenced for each population. Collapse the three columns and save the resulting dataset. And filter out the intergenic
data<- data %>% transmute(Gene = data$Gene, Details = data$Details, Population = data$Population, frequency= `F1 I1 R1` +`F1 I2 R1`+`F1 I3 R1`) %>% 
  filter(Details != "intergenic") %>%
  replace(is.na(.), 0)

singlegen_c_hyper("Deatherage2017", "glucose minimal medium", "Ecoli_K12")

### McCloskey2018
#trying to find a better way with a single function to do all of the following functions
#for gnd, filter out everything thats not gnd
data<- data %>% 
  filter(Population!="Evo04pgiEvo02EP",Population!="Evo04Evo01EP",Population!="Evo04Evo02EP",Population!="Evo04pgiEvo01EP",Population!="Evo04pgiEvo03EP",Population!="Evo04pgiEvo04EP",Population!="Evo04pgiEvo05EP",Population!="Evo04pgiEvo06EP",Population!="Evo04pgiEvo07EP",Population!="Evo04pgiEvo08EP",Population!="Evo04ptsHIcrrEvo01EP",Population!="Evo04ptsHIcrrEvo02EP",Population!="Evo04ptsHIcrrEvo03EP", Population!="Evo04ptsHIcrrEvo04EP",Population!= "Evo04sdhCBEvo02EP",Population!="Evo04sdhCBEvo03EP",Population!="Evo04tpiAEvo01EP",Population!="Evo04tpiAEvo02EP",Population!="Evo04tpiAEvo04EP",Population!="Evo04tpiAEvo03EP",Population!="Evo04sdhCBEvo01EP")

singlegen_c_hyper("McCloskey2018", "M9 minimal medium", "Ecoli_K12")

#for pgi, filter out everything thats not pgi
data<-data %>% 
  filter(Population!="Evo04Evo01EP",Population!="Evo04Evo02EP",Population!="Evo04gndEvo01EP",Population!="Evo04gndEvo02EP",Population!="Evo04gndEvo03EP",Population!="Evo04ptsHIcrrEvo01EP",Population!="Evo04ptsHIcrrEvo02EP",Population!="Evo04ptsHIcrrEvo03EP", Population!="Evo04ptsHIcrrEvo04EP",Population!= "Evo04sdhCBEvo02EP",Population!="Evo04sdhCBEvo03EP",Population!="Evo04tpiAEvo01EP",Population!="Evo04tpiAEvo02EP",Population!="Evo04tpiAEvo04EP",Population!="Evo04tpiAEvo03EP",Population!="Evo04sdhCBEvo01EP")

singlegen_c_hyper("McCloskey2018", "M9 minimal medium", "Ecoli_K12")

#for ptsHIcrr, filter out everything thats not ptsHIcrr
data<-data %>% 
  filter(Population!="Evo04pgiEvo02EP",Population!="Evo04Evo01EP",Population!="Evo04Evo02EP",Population!="Evo04gndEvo01EP",Population!="Evo04gndEvo02EP",Population!="Evo04gndEvo03EP",Population!="Evo04pgiEvo01EP",Population!="Evo04pgiEvo03EP",Population!="Evo04pgiEvo04EP",Population!="Evo04pgiEvo05EP",Population!="Evo04pgiEvo06EP",Population!="Evo04pgiEvo07EP",Population!="Evo04pgiEvo08EP",Population!= "Evo04sdhCBEvo02EP",Population!="Evo04sdhCBEvo03EP",Population!="Evo04tpiAEvo01EP",Population!="Evo04tpiAEvo02EP",Population!="Evo04tpiAEvo04EP",Population!="Evo04tpiAEvo03EP",Population!="Evo04sdhCBEvo01EP")

singlegen_c_hyper("McCloskey2018", "M9 minimal medium", "Ecoli_K12")

#for sdhCBE, filter out everything thats not sdhCBE
data<-data %>% 
  filter(Population!="Evo04pgiEvo02EP",Population!="Evo04Evo01EP",Population!="Evo04Evo02EP",Population!="Evo04gndEvo01EP",Population!="Evo04gndEvo02EP",Population!="Evo04gndEvo03EP",Population!="Evo04pgiEvo01EP",Population!="Evo04pgiEvo03EP",Population!="Evo04pgiEvo04EP",Population!="Evo04pgiEvo05EP",Population!="Evo04pgiEvo06EP",Population!="Evo04pgiEvo07EP",Population!="Evo04pgiEvo08EP",Population!="Evo04ptsHIcrrEvo01EP",Population!="Evo04ptsHIcrrEvo02EP",Population!="Evo04ptsHIcrrEvo03EP", Population!="Evo04ptsHIcrrEvo04EP",Population!="Evo04tpiAEvo01EP",Population!="Evo04tpiAEvo02EP",Population!="Evo04tpiAEvo04EP",Population!="Evo04tpiAEvo03EP")

singlegen_c_hyper("McCloskey2018", "M9 minimal medium", "Ecoli_K12")

#for tpiAE, filter out everything thats not tpiAE
data<-data %>% 
  filter(Population!="Evo04pgiEvo02EP",Population!="Evo04Evo01EP",Population!="Evo04Evo02EP",Population!="Evo04gndEvo01EP",Population!="Evo04gndEvo02EP",Population!="Evo04gndEvo03EP",Population!="Evo04pgiEvo01EP",Population!="Evo04pgiEvo03EP",Population!="Evo04pgiEvo04EP",Population!="Evo04pgiEvo05EP",Population!="Evo04pgiEvo06EP",Population!="Evo04pgiEvo07EP",Population!="Evo04pgiEvo08EP",Population!="Evo04ptsHIcrrEvo01EP",Population!="Evo04ptsHIcrrEvo02EP",Population!="Evo04ptsHIcrrEvo03EP", Population!="Evo04ptsHIcrrEvo04EP",Population!= "Evo04sdhCBEvo02EP",Population!="Evo04sdhCBEvo03EP",Population!="Evo04sdhCBEvo01EP")

singlegen_c_hyper("McCloskey2018", "M9 minimal medium", "Ecoli_K12")

#for the reference(evo), filter out everything thats not evo
data<-data %>% 
  filter(Population!="Evo04pgiEvo02EP",Population!="Evo04gndEvo01EP",Population!="Evo04gndEvo02EP",Population!="Evo04gndEvo03EP",Population!="Evo04pgiEvo01EP",Population!="Evo04pgiEvo03EP",Population!="Evo04pgiEvo04EP",Population!="Evo04pgiEvo05EP",Population!="Evo04pgiEvo06EP",Population!="Evo04pgiEvo07EP",Population!="Evo04pgiEvo08EP",Population!="Evo04ptsHIcrrEvo01EP",Population!="Evo04ptsHIcrrEvo02EP",Population!="Evo04ptsHIcrrEvo03EP", Population!="Evo04ptsHIcrrEvo04EP",Population!= "Evo04sdhCBEvo02EP",Population!="Evo04sdhCBEvo03EP",Population!="Evo04tpiAEvo01EP",Population!="Evo04tpiAEvo02EP",Population!="Evo04tpiAEvo04EP",Population!="Evo04tpiAEvo03EP",Population!="Evo04sdhCBEvo01EP")

singlegen_c_hyper("McCloskey2018", "M9 minimal medium", "Ecoli_K12")

#or just run as one function,each being an individiual population
singlegen_c_hyper("McCloskey2018", "M9 minimal medium", "Ecoli_K12")
