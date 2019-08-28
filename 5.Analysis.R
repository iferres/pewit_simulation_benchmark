PEWIT <- 'pewit_vs_sim/PEWIT.RDS'
ROARYI95 <- 'roaryi95_vs_sim/ROARY_I95.RDS'
ROARYI70 <- 'roaryi70_vs_sim/ROARY_I70.RDS'
PANX <- 'panx_vs_sim/PANX.RDS'
MICROPANBLASTSINGLE <- 'micropanBlast_vs_sim/MICROPAN_BLAST_SINGLE.RDS'
MICROPANBLASTCOMPLETE <- 'micropanBlast_vs_sim/MICROPAN_BLAST_COMPLETE.RDS'
MICROPANBLASTAVERAGE <- 'micropanBlast_vs_sim/MICROPAN_BLAST_AVERAGE.RDS'

rds <- c(PEWIT = PEWIT,
         ROARYI70 = ROARYI95, 
         ROARYI95 = ROARYI70,
         PANX = PANX,
         MICROPANBLASTSINGLE = MICROPANBLASTSINGLE,
         MICROPANBLASTCOMPLETE = MICROPANBLASTCOMPLETE,
         MICROPANBLASTAVERAGE = MICROPANBLASTAVERAGE)


res <- lapply(rds, readRDS)

res <- do.call(rbind, res)
res <- res[-which(rowSums(res[, 1:4])==0), ]

################
## Statistics ##
################


# Sensitivity, or True Positive Rate, or Recall
res$Recall <- res$TP / (res$TP + res$FN)

# Specificity, or Precision
res$Precision <- res$TP / (res$TP + res$FP)

# False Positive Rate
res$FPR <- res$FP / (res$FP + res$TN)

# F1 Score
res$F1_Score <- 2 * ((res$Precision * res$Recall) / 
                       (res$Precision + res$Recall))
# 0 / 0 == NaN; but:
# (x*y) / (x+y) --> 0 when x and y --> 0
# .. so:
res$F1_Score[is.nan(res$F1_Score)] <- 0



##########
## Plot ##
##########

res$Ne_F <- as.factor(res$Ne)


library(ggpubr)

# ROC curve
ggscatter(res, 
          x='FPR', 
          y='Recall', 
          color='Software', 
          facet.by = 'Ne',
          # shape='Ne_F', 
          ellipse = T, 
          ellipse.type = 'convex')


# F1 Score
ggboxplot(res, 
          x = 'Software', 
          y = 'F1_Score', 
          facet.by = 'Ne_F', 
          scales='free_y', 
          color = 'Software')

ggline(res, 
       x = 'Ne', 
       y = 'F1_Score',  
       scales='free_y', add = c('median'),  
       color = 'Software')


library(ggthemes)
library(ggsci)
ggboxplot(res, x = 'Software',
       y = 'F1_Score',
       scales='free_y', 
       # add = c('boxplot'),  
       facet.by = 'Ne',
       color = 'Software')  + 
  theme_tufte() + 
  ggsci::scale_color_npg() + 
  rremove('x.text') + 
  rremove('x.ticks')

