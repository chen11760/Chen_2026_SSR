##############################################################################################################
# Replicating figures for project:
# Reassessing the Gendered Link Between Maternal Employment and Adult Children's Labor Market Participations: 
# A Trajectory-based Approach at Social Science Research

# Data cleaning steps for NLSY79 & CNLSY79 not included

# Updated by Xueqian(Chelsea) Chen
# 18th June 2026
##############################################################################################################

packages <- c("haven","dplyr","RColorBrewer","TraMineR","gridExtra","cluster",
              "foreign","factoextra","ggplot2","WeightedCluster","tidyr","readxl",
              "lme4","plm","clubSandwich","broom.mixed","lmerTest","margins",
              "ggeffects","cowplot","magick")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

setwd("/Users/xueqian/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Research/maternal employment intergenerational/data_sa_firstbirth")


##############################################################################################################
# Sequence analysis and clustering
##############################################################################################################

### Import cleaned data for employment sequence
motherwork <- read_dta("motherwork1.dta")

### Define categories for work status
tracwork <- motherwork %>% 
  select(starts_with('work')) 

labs <- c("Full-time Employed","Part-time Employed","Marginally Employed","On Leave With Work", "Employed With Unknown Time", "Not Employed")

palette <- brewer.pal(length(labs), 'RdYlBu')

seqwork <- seqdef(tracwork, labels=labs, cpal=palette)


### Caluclate number of distinct sequences: 
seqtab(seqwork, idx=0) %>% nrow  # n = 1,889


### Calculate pairwise distance matrix using Dynamic Hamming distance
dis <- seqdist(seqwork, method = "DHD", norm = "auto", sm=NULL)



### Determining typology of sequences

# Hierarchical agglomerative cluster
agnes <- as.dist(dis) %>% agnes(method="ward", keep.diss=FALSE)

# dendrogram plot
plot(as.dendrogram(agnes),leaflab="none")
abline(h = 5, col = 'red')

# Indicators of quality
wardRange <- as.clustrange(agnes, diss=dis, ncluster=6)
wardRange 
summary(wardRange, max.rank=2)  # A 3-cluster solution being the best, followed by a 5-cluster solution

# Partition by 3 or 5
part <- cutree(agnes, 3) 
partnew <- cutree(agnes,5) 

# Number of members in each solution
table(part)
table(partnew)



### Plot sequence by type for each solution

# 3-cluster solution
prop <-round(prop.table(table(part))*100, digits = 1)
par(mfrow=c(2,2), mar=c(2.5,1.8,1.8,1.8))
for(i in 1:3) seqdplot(seqwork[part==i,], xtlab=1:216,border=NA,with.legend=FALSE, main=paste('cluster',i,'(',prop[i],'%)',sep=''))
seqlegend(seqwork, cex=.8)  
# 3-cluster solution: mostly not employed, mostly full-time employed, mostly part-time employed


# 5-cluster solution
propnew <-round(prop.table(table(partnew))*100, digits = 1)
par(mfrow=c(2,3), mar=c(2.5,1.8,1.8,1.8))
for(i in 1:5) seqdplot(seqwork[partnew==i,], xtlab=1:216,border=NA,with.legend=FALSE, main=paste('cluster',i,'(',propnew[i],'%)',sep=''))
seqlegend(seqwork, cex=.8)
# 5-cluster solution include: late return to work, early return to work, mostly full-time employed, mostly part-time employed, mostly not employed
# for the purpose of this study (to explore more nuances in maternal employment timing and continuity),
# I opted into 5-cluster solution.






##############################################################################################################
# FIGURE 1. PATTERNS OF MATERNAL EMPLOYMENT TRAJECTORIES (NLSY79)
##############################################################################################################

### Figure1a: index plot
png("Figure1a.png",
     width = 4000,
     height = 3000,
     res = 600)

trajlabs <- c("Late Return","Early Return","Full-time","Part-time","Not Working")

par(mfrow=c(2,3), mar=c(2.5,1.8,1.8,1.8))

seqIplot(seqwork[partnew==3,], xtlab=FALSE,border=NA,with.legend=FALSE, main=paste(1,":", trajlabs[3],'(',propnew[3],'%)',sep=''))
axis(1,at = c(1,37,73,109,145,181,216), labels = c("0","3","6","9","12","15","18"))

seqIplot(seqwork[partnew==2,], xtlab=FALSE,border=NA,with.legend=FALSE, main=paste(2,":", trajlabs[2],'(',propnew[2],'%)',sep=''))
axis(1,at = c(1,37,73,109,145,181,216), labels = c("0","3","6","9","12","15","18"))

seqIplot(seqwork[partnew==1,], xtlab=FALSE,border=NA,with.legend=FALSE, main=paste(3,":", trajlabs[1],'(',propnew[1],'%)',sep=''))
axis(1,at = c(1,37,73,109,145,181,216), labels = c("0","3","6","9","12","15","18"))

seqIplot(seqwork[partnew==4,], xtlab=FALSE,border=NA,with.legend=FALSE, main=paste(4,":", trajlabs[4],'(',propnew[4],'%)',sep=''))
axis(1,at = c(1,37,73,109,145,181,216), labels = c("0","3","6","9","12","15","18"))

seqIplot(seqwork[partnew==5,], xtlab=FALSE,border=NA,with.legend=FALSE, main=paste(5,":", trajlabs[5],'(',propnew[5],'%)',sep=''))
axis(1,at = c(1,37,73,109,145,181,216), labels = c("0","3","6","9","12","15","18"))

seqlegend(seqwork, cex=0.8)

dev.off()




### Figure1b: distribution plot
png(
  "Figure1b.png",
  width = 4000,
  height = 3000,
  res = 600
)

trajlabs <- c("Late Return","Early Return","Full-time","Part-time","Not Working")

par(mfrow=c(2,3), mar=c(2.5,1.8,1.8,1.8))

seqdplot(seqwork[partnew==3,], xtlab=FALSE,border=NA, with.legend=FALSE, main=paste(1,":", trajlabs[3],'(',propnew[3],'%)',sep=''))
axis(1,at = c(0,36,73,109,145,181,216), labels = c("0","3","6","9","12","15","18"))

seqdplot(seqwork[partnew==2,], xtlab=FALSE, border=NA, with.legend=FALSE, main=paste(2,":", trajlabs[2],'(',propnew[2],'%)',sep=''))
axis(1,at = c(1,37,73,109,145,181,216), labels = c("0","3","6","9","12","15","18"))

seqdplot(seqwork[partnew==1,], xtlab=FALSE, border=NA, with.legend=FALSE, main=paste(3,":", trajlabs[1],'(',propnew[1],'%)',sep=''))
axis(1,at = c(1,37,73,109,145,181,216), labels = c("0","3","6","9","12","15","18"))

seqdplot(seqwork[partnew==4,], xtlab=FALSE, border=NA, with.legend=FALSE, main=paste(4,":", trajlabs[4],'(',propnew[4],'%)',sep=''))
axis(1,at = c(1,37,73,109,145,181,216), labels = c("0","3","6","9","12","15","18"))

seqdplot(seqwork[partnew==5,], xtlab=FALSE, border=NA, with.legend=FALSE, main=paste(5,":", trajlabs[5],'(',propnew[5],'%)',sep=''))
axis(1,at = c(1,37,73,109,145,181,216), labels = c("0","3","6","9","12","15","18"))

seqlegend(seqwork, cex=0.8)

dev.off()





##############################################################################################################
# FIGURE 2. DISTRIBUTION OF MATERNAL EMPLOYMENT STATES BY TRAJECTORY PATTERNS
##############################################################################################################

# Import summarized data for Figure 2
figure2 <- read_excel("des1.xlsx")

# Transform to long data
figure2_long <- figure2 %>%
  pivot_longer(cols = starts_with("n_"), 
               names_to = "Work_Status", 
               values_to = "Average_Months")
figure2_long$Work_Status <- factor(figure2_long$Work_Status, 
                              levels = c("n_full", "n_part", "n_marginal", "n_leave", "n_unknown", "n_no"),
                              labels = c("Full-time", "Part-time", "Marginally", "On Leave", "Unknown", "Not Working"))

# reorder
figure2_long$Work_Trajectory <- factor(
  figure2_long$Work_Trajectory,
  levels = rev(c("Full-time",
                 "Early Return",
                 "Late Return",
                 "Part-time",
                 "Not Working")))

# Plot
palette <- brewer.pal(6, "RdYlBu") 
ggplot(figure2_long, aes(x = Work_Trajectory, y = Average_Months, fill = Work_Status)) +
  geom_bar(stat = "identity", position = "stack") +  
      labs(x = "",
           y = "Average Months",
           fill = "Work Status") +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 12, face = "bold"),
        axis.text.y  = element_text(size = 12, face = "bold")) +
  scale_fill_manual(values = c("Full-time" = palette[1], 
                               "Part-time" = palette[2], 
                               "Marginally" = palette[3], 
                               "On Leave" = palette[4], 
                               "Unknown" = palette[5], 
                               "Not Working" = palette[6]) ) +
  coord_flip()

ggsave("Figure2.png", width = 10, height = 5, dpi = 600)
