# Roary's publication ref file
# Kallonen, Teemu (2017): Pan genome reference fasta. figshare. Dataset.
# https://doi.org/10.6084/m9.figshare.4873022.v1
# 'https://figshare.com/articles/Pan_genome_reference_fasta/4873022'
ref <- 'Ecoli_pan_genome_reference.fa'
if (!file.exists(ref)){
  stop('You have to manually downlowad pangenome reference file, name it 
       "Ecoli_pan_genome_reference.fa", and place it on working directory.')
}

library(simurg)

# Simulated evolution times
nes <- c(2e+09, 1e+10, 5e+10, 25e+10, 125e+10)

# Replicates
reps <- 5

set.seed(123)
for (n in nes){
  
  for (r in seq_len(reps)){
    
    Nexp <- sub('e[+]', 'E', n)
    id <- paste0(Nexp, '_rep', r)
    dout <- paste('sim', id, sep = '_')
    pg <- simpg(ref = ref, 
                norg = 10, 
                ne = n, 
                C = 100, 
                u = 1E-8, 
                v = 1E-11, 
                mu = 5E-12, 
                write_by = 'genome', 
                dir_out = dout, 
                replace = TRUE, 
                verbose = FALSE)
    rds <- paste('pg_', id, '.RDS', sep = '')
    saveRDS(pg, file = paste(dout, rds, sep = '/'))
    message(rds)
    
  }
  
}

## NEXT: format_input.R
