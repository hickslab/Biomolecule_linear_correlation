library(devtools)
library(tidyverse)

setwd("Z:/Lab_Members/Patric/Projects/Correlation/Diss/Markdown_SOP/")


##### read in the data #####
mdata <- read_csv("metabolite_data.csv")
pdata <- read_csv("protein_data.csv")


#format the data
metsum <- mdata %>% 
  gather(key = "Sample", value = "MAUC", 2:length(mdata)) %>% 
  select(Mass, Sample, MAUC)

#format the proteomics data

prosum <- pdata %>% 
  gather(key = "Sample", value = "PAUC", 2:length(pdata)) %>% 
  select(Accession, Sample, PAUC) 

MPdf <- full_join(metsum, prosum)

MPdf$Mass <- as.factor(MPdf$Mass)

MPdf <- MPdf %>% 
  select(Mass, Sample, Accession, MAUC, PAUC) %>% 
  separate(Sample, into = c("Sample", "Replicate"), sep = "-")


##### Define the function to loop through metabolites #####


cor_prot_metabolites <- function(MPdf, output){
  
  proteins <- unique(MPdf$Accession) #get a list of each protein in the dataframe
  metabolites <- unique(MPdf$Mass) #get a list of the metabolites in the dataframe
  output <- data.frame()
  
  for(j in 1:length(metabolites)){
    
    temp1 <- MPdf %>% filter(Mass == metabolites[j])
    
    for(i in 1:length(proteins)){
      temp <- temp1 %>% 
        filter(Accession == proteins[i])
      
      temp.pearson <- cor.test(temp$MAUC, temp$PAUC, method = "pearson")                 
      
      temp$pearson.cor <- temp.pearson$estimate
      temp$pearson.sig <- temp.pearson$p.value
      
      temp.spearman <- cor.test(temp$MAUC, temp$PAUC, method = "spearman", exact = FALSE)                 
      
      temp$spearman.rho <- temp.spearman$estimate
      temp$spearman.sig <- temp.spearman$p.value
      
      output <- rbind.data.frame(output, temp)
    }}
  return(output)
}

MPdf_corr <- cor_prot_metabolites(MPdf = MPdf)

#reduce the data to fewer numbers
MPdf_corr2 <-  MPdf_corr %>% 
  group_by(Mass, Sample, Accession) %>% 
  summarize("mMAUC" = mean(MAUC),
            "sdMAUC"= sd(MAUC),
            "mPAUC" = mean(PAUC),
            "sdPAUC"= sd(PAUC),
            "pearson.cor" = mean(pearson.cor),
            "pearson.sig" = mean(pearson.sig),
            "spearman.rho" = mean(spearman.rho),
            "spearman.sig" = mean(spearman.sig))

MPdf_corr3 <- MPdf_corr2 %>% 
  filter(Sample == "ConditionA")

#plots

theme_custom <- function(base_size = 8){
  theme_bw(base_size = base_size) %+replace%
    theme(
      strip.background = element_blank(),
      #axis.ticks =  element_line(colour = "black"),
      panel.background = element_blank(),
      #panel.border = element_blank(),
      #panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      plot.background = element_blank(),
      #plot.margin = unit(c(0.5,  0.5, 0.5, 0.5), "lines"),
      #axis.line.x = element_line(color="black", size = 1),
      #axis.line.y = element_line(color="black", size = 1)
    )
}




#graph a jitter of the pearson and spearman without significance
MPdf_corr3 %>% 
  ungroup() %>% 
  select(Mass, Accession, pearson.cor, spearman.rho, pearson.sig, spearman.sig) %>% 
  gather(key = "test", value = "sig", 3:6) %>% 
  filter(test == "pearson.cor" | test == "spearman.rho") %>% 
  ggplot(., aes(x = test, y = sig, color = Mass))+
  geom_jitter(width = 0.1, alpha = 0.8, size = 0.5)+
  geom_boxplot(fill = NA, outlier.shape = NA)+
  facet_wrap(.~Mass, nrow = 2)+
  theme_custom()+
  scale_y_continuous(limits = c(-1, 1))+
  coord_flip()+
  scale_color_viridis_d(begin = 0.3, end = 0.7)+
  labs(x = "Coefficient",
       y = "Correlation test")

