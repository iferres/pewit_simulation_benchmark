# setwd("/export/home/iferres/Escritorio/simulations_pewit_eval")

dins <- dir(pattern = 'sim_.+_rep\\d+')
names(dins) <- dins
source('3.Benchmark-Functions.R')

# Sounds good, doesn't work:
# library(parallel)
# res_roary_i70 <- mclapply(dins, bench_roary_i70, mc.cores = 2)

#################
### ROARY I70 ###
#################

res <- data.frame(TN = integer(),
                  FP = integer(),
                  FN = integer(),
                  TP = integer(),
                  Software = character(),
                  Ne = numeric(),
                  Replicate = integer(), 
                  stringsAsFactors = FALSE)

for (i in seq_along(dins)){
  spl <- strsplit(dins[i], '_')[[1]]
  res[i, c('TN','FP','FN','TP')] <- unclass(bench_roary_i70(dins[i]))
  res[i, 'Software'] <- 'roary_i70'
  res[i, 'Ne'] <- as.numeric(spl[2])
  res[i, 'Replicate'] <- as.integer(sub('rep', '', spl[3]))
}

saveRDS(res, file = 'roaryi70_vs_sim/ROARY_I70.RDS')


#################
### ROARY I95 ###
#################

res <- data.frame(TN = integer(),
                  FP = integer(),
                  FN = integer(),
                  TP = integer(),
                  Software = character(),
                  Ne = numeric(),
                  Replicate = integer(), 
                  stringsAsFactors = FALSE)

for (i in seq_along(dins)){
  spl <- strsplit(dins[i], '_')[[1]]
  res[i, c('TN','FP','FN','TP')] <- unclass(bench_roary_i95(dins[i]))
  res[i, 'Software'] <- 'roary_i95'
  res[i, 'Ne'] <- as.numeric(spl[2])
  res[i, 'Replicate'] <- as.integer(sub('rep', '', spl[3]))
}

saveRDS(res, file = 'roaryi95_vs_sim/ROARY_I95.RDS')


############
### PANX ###
############

res <- data.frame(TN = integer(),
                  FP = integer(),
                  FN = integer(),
                  TP = integer(),
                  Software = character(),
                  Ne = numeric(),
                  Replicate = integer(), 
                  stringsAsFactors = FALSE)

for (i in seq_along(dins)){
  spl <- strsplit(dins[i], '_')[[1]]
  res[i, c('TN','FP','FN','TP')] <- unclass(bench_panx(dins[i]))
  res[i, 'Software'] <- 'PanX'
  res[i, 'Ne'] <- as.numeric(spl[2])
  res[i, 'Replicate'] <- as.integer(sub('rep', '', spl[3]))
}

saveRDS(res, file = 'panx_vs_sim/PANX.RDS')




#############
### PEWIT ###
#############

res <- data.frame(TN = integer(),
                  FP = integer(),
                  FN = integer(),
                  TP = integer(),
                  Software = character(),
                  Ne = numeric(),
                  Replicate = integer(), 
                  stringsAsFactors = FALSE)

for (i in seq_along(dins)){
  spl <- strsplit(dins[i], '_')[[1]]
  res[i, c('TN','FP','FN','TP')] <- unclass(bench_pewit(dins[i], 
                                                        hmm_pfam = '~/Escritorio/pfam-A/Pfam-A.hmm',
                                                        dat_pfam = '~/Escritorio/pfam-A/Pfam-A.hmm.dat'))
  res[i, 'Software'] <- 'Pewit'
  res[i, 'Ne'] <- as.numeric(spl[2])
  res[i, 'Replicate'] <- as.integer(sub('rep', '', spl[3]))
}

saveRDS(res, file = 'pewit_vs_sim/PEWIT.RDS')



######################
### MICROPAN BLAST ###
######################

# Single Linkage (micropan default)

res <- data.frame(TN = integer(),
                  FP = integer(),
                  FN = integer(),
                  TP = integer(),
                  Software = character(),
                  Ne = numeric(),
                  Replicate = integer(), 
                  stringsAsFactors = FALSE)

for (i in seq_along(dins)){
  spl <- strsplit(dins[i], '_')[[1]]
  res[i, c('TN','FP','FN','TP')] <- unclass(bench_micropan_blast(dins[i]))
  res[i, 'Software'] <- 'Micropan_BlastAllAll_Single'
  res[i, 'Ne'] <- as.numeric(spl[2])
  res[i, 'Replicate'] <- as.integer(sub('rep', '', spl[3]))
}

saveRDS(res, file = 'micropanBlast_vs_sim/MICROPAN_BLAST_SINGLE.RDS')

# The following asumes that single linkage micropan has already been ran,
# it uses the already computed blastAllAll result.

# Complete Linkage

res <- data.frame(TN = integer(),
                  FP = integer(),
                  FN = integer(),
                  TP = integer(),
                  Software = character(),
                  Ne = numeric(),
                  Replicate = integer(), 
                  stringsAsFactors = FALSE)

for (i in seq_along(dins)){
  spl <- strsplit(dins[i], '_')[[1]]
  res[i, c('TN','FP','FN','TP')] <- unclass(bench_micropan_blast_complete(dins[i]))
  res[i, 'Software'] <- 'Micropan_BlastAllAll_Complete'
  res[i, 'Ne'] <- as.numeric(spl[2])
  res[i, 'Replicate'] <- as.integer(sub('rep', '', spl[3]))
}

saveRDS(res, file = 'micropanBlast_vs_sim/MICROPAN_BLAST_COMPLETE.RDS')

# Average Linkage

res <- data.frame(TN = integer(),
                  FP = integer(),
                  FN = integer(),
                  TP = integer(),
                  Software = character(),
                  Ne = numeric(),
                  Replicate = integer(), 
                  stringsAsFactors = FALSE)

for (i in seq_along(dins)){
  spl <- strsplit(dins[i], '_')[[1]]
  res[i, c('TN','FP','FN','TP')] <- unclass(bench_micropan_blast_average(dins[i]))
  res[i, 'Software'] <- 'Micropan_BlastAllAll_Average'
  res[i, 'Ne'] <- as.numeric(spl[2])
  res[i, 'Replicate'] <- as.integer(sub('rep', '', spl[3]))
}

saveRDS(res, file = 'micropanBlast_vs_sim/MICROPAN_BLAST_AVERAGE.RDS')



####################
## MICROPAN HMMER ##
####################

res <- data.frame(TN = integer(),
                  FP = integer(),
                  FN = integer(),
                  TP = integer(),
                  Software = character(),
                  Ne = numeric(),
                  Replicate = integer(), 
                  stringsAsFactors = FALSE)

for (i in seq_along(dins)){
  spl <- strsplit(dins[i], '_')[[1]]
  res[i, c('TN','FP','FN','TP')] <- unclass(bench_micropan_hmmer(dins[i], 
                                                                 hmm_pfam = '~/Escritorio/pfam-A/Pfam-A.hmm'))
  res[i, 'Software'] <- 'Micropan_Hmmer'
  res[i, 'Ne'] <- as.numeric(spl[2])
  res[i, 'Replicate'] <- as.integer(sub('rep', '', spl[3]))
}

saveRDS(res, file = 'micropanHMMER_vs_sim/MICROPAN_HMMER.RDS')

