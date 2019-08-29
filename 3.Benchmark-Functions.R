

# Compute contingency table function
`%contingency%` <- function(x, y){
  
  #Order
  dn <- dimnames(x)[[1]]
  x <- x[dn, dn]
  y <- y[dn, dn]
  
  #Clean
  diag(x) <- 0L
  diag(y) <- 0L
  
  #Compute contigency
  # tbl <- table(x%checkv%y)
  tbl <- table(factor(paste(x, y, sep = ''), levels = c('00', '01', '10', '11')))
  names(tbl)[which(names(tbl)==c('00', '01', '10', '11'))] <- c('TN','FP','FN','TP')
  #Substract diagonal of "TRUE NEGATIVES"
  dm <- dim(x)
  tbl['TN'] <- tbl['TN'] - dm[2]
  #Values are repeated, so..
  res <- tbl/2
  
  return(res)
}


library(pagoo)

bench_roary_i70 <- function(din, dout = 'roaryi70_vs_sim', cmd = 'singularity run singularity_containers/roary.sif'){
  
  gffs <- list.files(path = din, pattern = '[.]gff$', full.names = TRUE)
  sim <- list.files(path = din, pattern = '[.]RDS$', full.names = TRUE)
  
  if (!dir.exists(dout)){
    dir.create(dout)
  }
  # dn <- dirname(gffs[1])
  f <- paste0(dout, '/', din)
  if (dir.exists(f)) unlink(f, recursive = TRUE)
  run <- paste(cmd, '-p 4 -i 70 -f', f, paste(gffs, collapse = ' '), '2> /dev/null')
  system(run)
  
  # Process output
  # Load simulation
  simu <- readRDS(sim)
  clus <- simu$gene_list
  seqs <- unlist(clus)
  nseqs <- length(seqs)
  mat_expect <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
  for (i in 1:length(clus)){
    ge <- clus[[i]]
    mat_expect[ge, ge] <- 1L
  }
  
  
  # Load run
  csv <- paste(f, '/gene_presence_absence.csv', sep = '')
  pg <- pagoo::roary_2_pagoo(csv)
  grp <- split(as.character(pg$.__enclos_env__$private$.DF$gene), 
               pg$.__enclos_env__$private$.DF$group)
  mat_observed <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
  for (i in 1:length(grp)){
    ge <- grp[[i]]
    mat_observed[ge, ge] <- 1L
  }
  
  
  # Compute Contingency
  res <- mat_expect %contingency% mat_observed
  return(res)
}



bench_roary_i95 <- function(din, dout = 'roaryi95_vs_sim', cmd = 'singularity run singularity_containers/roary.sif'){
  
  gffs <- list.files(path = din, pattern = '[.]gff$', full.names = TRUE)
  sim <- list.files(path = din, pattern = '[.]RDS$', full.names = TRUE)
  
  if (!dir.exists(dout)){
    dir.create(dout)
  }
  # dn <- dirname(gffs[1])
  f <- paste0(dout, '/', din)
  if (dir.exists(f)) unlink(f, recursive = TRUE)
  run <- paste(cmd, '-p 4 -i 95 -f', f, paste(gffs, collapse = ' '), '2> /dev/null')
  system(run)
  
  # Process output
  # Load simulation
  simu <- readRDS(sim)
  clus <- simu$gene_list
  seqs <- unlist(clus)
  nseqs <- length(seqs)
  mat_expect <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
  for (i in 1:length(clus)){
    ge <- clus[[i]]
    mat_expect[ge, ge] <- 1L
  }
  
  
  # Load run
  csv <- paste(f, '/gene_presence_absence.csv', sep = '')
  pg <- pagoo::roary_2_pagoo(csv)
  grp <- split(as.character(pg$.__enclos_env__$private$.DF$gene), 
               pg$.__enclos_env__$private$.DF$group)
  mat_observed <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
  for (i in 1:length(grp)){
    ge <- grp[[i]]
    mat_observed[ge, ge] <- 1L
  }
  
  
  # Compute Contingency
  res <- mat_expect %contingency% mat_observed
  return(res)
}



bench_panx <- function(din, dout = 'panx_vs_sim', cmd = 'singularity run singularity_containers/panx.sif'){
  gbks <- list.files(path = din, pattern = '[.]gbk$', full.names = TRUE)
  sim <- list.files(path = din, pattern = '[.]RDS$', full.names = TRUE)
  if (!dir.exists(dout)){
    dir.create(dout)
  }
  
  fn <- paste0(dout, '/', din)
  if (dir.exists(fn)){
    unlink(fn, recursive = TRUE)
  }
  dir.create(fn)
  file.symlink(normalizePath(gbks), paste(fn, basename(gbks), sep = '/'))
  
  run <- paste(cmd, '-fn', fn, '-sl Simulated_bacteria -cg 0.7 -st 1 2 3 4 5 6 -t 4 -ct')
  tr <- try(system(run))
  
  #There may have problems parsing gbk files
  if (!tr){
    
    #clean: remove all files but allclusters_final.tsv and allclusters.tsv
    alls <- list.files(path = fn, recursive = T, full.names = TRUE)
    torm <- grep('allclusters.*[.]tsv$', alls, invert = TRUE, value = TRUE)
    file.remove(torm)
    
    tsvs <- dir(path = fn, recursive = TRUE)
    tsv2 <- paste(fn, basename(tsvs), sep = '/')
    file.rename(paste(fn, tsvs, sep = '/'), tsv2)
    unlink(list.dirs(fn)[-1], recursive = T)
    gpl <- grepl('allclusters_final.tsv', tsv2, fixed = TRUE)
    if (any(gpl)){
      allclusters_final_tsv <- tsv2[which(gpl)]
    }else{
      # stop('No allclusters_final.tsv file.')
      allclusters_final_tsv <- tsv2[1]
    }
    
    # Process output
    # Load simulation
    simu <- readRDS(sim)
    clus <- simu$gene_list
    seqs <- unlist(clus)
    nseqs <- length(seqs)
    mat_expect <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
    for (i in 1:length(clus)){
      ge <- clus[[i]]
      mat_expect[ge, ge] <- 1L
    }
    
    pg <- pagoo:::panx_2_pagoo(allclusters_final_tsv = allclusters_final_tsv, sep = '_')
    grp <- split(as.character(pg$.__enclos_env__$private$.DF$gene), 
                 pg$.__enclos_env__$private$.DF$group)
    mat_observed <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
    for (i in 1:length(grp)){
      ge <- grp[[i]]
      mat_observed[ge, ge] <- 1L
    }
    
    
    # Compute Contingency
    res <- mat_expect %contingency% mat_observed
    return(res)
    
    
  }else{
    return(as.table(c(TN = 0L, FP = 0L, FN = 0L, TP = 0L)))
  }
  
}





bench_pewit <- function(din, dout = 'pewit_vs_sim', hmm_pfam, dat_pfam){
  
  require(pewit)
  
  gffs <- list.files(path = din, pattern = '[.]gff$', full.names = TRUE)
  sim <- list.files(path = din, pattern = '[.]RDS$', full.names = TRUE)
  
  if (!dir.exists(dout)){
    dir.create(dout)
  }
  # dn <- dirname(gffs[1])
  dd <- paste0(dout, '/', din)
  if (dir.exists(dd)) unlink(dd, recursive = TRUE)
  dir.create(dd)
  
  pg <- pangenome(gffs = gffs, 
                  hmm_pfam = hmm_pfam, 
                  dat_pfam = dat_pfam, 
                  n_threads = 4)
  saveRDS(pg, file = paste(dd, 'pewit_out.RDS', sep = '/'))
  
  # Process output
  # Load simulation
  simu <- readRDS(sim)
  clus <- simu$gene_list
  seqs <- unlist(clus)
  nseqs <- length(seqs)
  mat_expect <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
  for (i in 1:length(clus)){
    ge <- clus[[i]]
    mat_expect[ge, ge] <- 1L
  }
  
  grp <- split(as.character(pg$.__enclos_env__$private$.DF$gene), 
               pg$.__enclos_env__$private$.DF$group)
  mat_observed <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
  for (i in 1:length(grp)){
    ge <- grp[[i]]
    mat_observed[ge, ge] <- 1L
  }
  
  
  # Compute Contingency
  res <- mat_expect %contingency% mat_observed
  return(res)
  
}



bench_micropan_blast <- function(din, dout = 'micropanBlast_vs_sim'){
  
  require(micropan)
  require(seqinr)
  
  faas <- list.files(path = din, pattern = '[.]faa$', full.names = TRUE)
  sim <- list.files(path = din, pattern = '[.]RDS$', full.names = TRUE)
  
  if (!dir.exists(dout)){
    dir.create(dout)
  }
  # dn <- dirname(gffs[1])
  dd <- paste0(dout, '/', din)
  if (dir.exists(dd)) unlink(dd, recursive = TRUE)
  dir.create(dd)
  
  
  gid <- paste0("%0",nchar(length(faas)),"d")
  gid <- paste0('GID',sprintf(gid,1:length(faas)))
  ref <- cbind(basename(faas), gid)
  write.csv(ref, file = paste0(dd, '/ref_gid.csv'),quote = F)
  
  prep_faas <- sapply(faas, function(x){
    ff <- basename(x)
    gd <- ref[which(ref[, 1] == ff), 2]
    out <- paste(dd, ff, sep = '/')
    panPrep(in.file = x, GID.tag = gd, out.file = out)
  })
  
  nams <- lapply(prep_faas, function(x){
    fa <- read.fasta(x)
    an <- sapply(fa, attr, 'Annot')
    spl <- strsplit(sub('^>', '', an), ' ')
    mm <- do.call(rbind, spl)
    rownames(mm) <- NULL
    mm
  })
  
  nams <- do.call(rbind, nams)
  
  blout_dir <- paste(dd, 'blast_out', sep = '/')
  if (dir.exists(blout_dir)) unlink(blout_dir, recursive = TRUE)
  dir.create(blout_dir)
  
  runMicropanBlast <- function(prep_faas, oblf){
    njb <- paste0(sprintf('%02d',sample(1:20, 2)), collapse = '')
    blastAllAll(prep_faas, 
                oblf, 
                e.value = 1e-10,
                threads = 4L,
                job = as.integer(njb),
                verbose = FALSE)
    bl <- list.files(path = oblf, full.names = TRUE)
    df <- bDist(bl, verbose = FALSE)
    # unlink(oblf, recursive = TRUE)
    bcl <- bClust(df, threshold = 0.75)
    return(bcl)
  }
  
  clust <- runMicropanBlast(prep_faas = prep_faas, oblf = blout_dir)
  names(clust) <- sapply(names(clust), function(x) nams[which(nams[,1]==x), 2])
  grp <- split(names(clust), clust)
  
  simu <- readRDS(sim)
  clus <- simu$gene_list
  seqs <- unlist(clus)
  nseqs <- length(seqs)
  mat_expect <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
  for (i in 1:length(clus)){
    ge <- clus[[i]]
    mat_expect[ge, ge] <- 1L
  }
  
  mat_observed <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
  for (i in 1:length(grp)){
    ge <- grp[[i]]
    mat_observed[ge, ge] <- 1L
  }
  
  
  res <- mat_expect %contingency% mat_observed
  return(res)
  
}

# bench_micropan_blast_single <- function(din, dout = 'micropanBlast_vs_sim'){
#   sim <- list.files(path = din, pattern = '[.]RDS$', full.names = TRUE)
#   
#   dd <- paste0(dout, '/', din)
#   prep_faas <- list.files(dd, pattern = '[.]faa$', full.names = TRUE)
#   nams <- lapply(prep_faas, function(x){
#     fa <- read.fasta(x)
#     an <- sapply(fa, attr, 'Annot')
#     spl <- strsplit(sub('^>', '', an), ' ')
#     mm <- do.call(rbind, spl)
#     rownames(mm) <- NULL
#     mm
#   })
#   
#   nams <- do.call(rbind, nams)
#   
#   blout_dir <- paste(dd, 'blast_out', sep = '/')
#   bl <- list.files(path = blout_dir, full.names = TRUE)
#   df <- bDist(bl, verbose = FALSE)
#   clust <- bClust(df, linkage = 'single', threshold = 0.75)
#   
#   names(clust) <- sapply(names(clust), function(x) nams[which(nams[,1]==x), 2])
#   grp <- split(names(clust), clust)
#   
#   simu <- readRDS(sim)
#   clus <- simu$gene_list
#   seqs <- unlist(clus)
#   nseqs <- length(seqs)
#   mat_expect <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
#   for (i in 1:length(clus)){
#     ge <- clus[[i]]
#     mat_expect[ge, ge] <- 1L
#   }
#   
#   mat_observed <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
#   for (i in 1:length(grp)){
#     ge <- grp[[i]]
#     mat_observed[ge, ge] <- 1L
#   }
#   
#   
#   res <- mat_expect %contingency% mat_observed
#   return(res)
# }


bench_micropan_blast_complete <- function(din, dout = 'micropanBlast_vs_sim'){
  sim <- list.files(path = din, pattern = '[.]RDS$', full.names = TRUE)
  
  dd <- paste0(dout, '/', din)
  prep_faas <- list.files(dd, pattern = '[.]faa$', full.names = TRUE)
  nams <- lapply(prep_faas, function(x){
    fa <- read.fasta(x)
    an <- sapply(fa, attr, 'Annot')
    spl <- strsplit(sub('^>', '', an), ' ')
    mm <- do.call(rbind, spl)
    rownames(mm) <- NULL
    mm
  })
  
  nams <- do.call(rbind, nams)
  
  blout_dir <- paste(dd, 'blast_out', sep = '/')
  bl <- list.files(path = blout_dir, full.names = TRUE)
  df <- bDist(bl, verbose = FALSE)
  clust <- bClust(df, linkage = 'complete', threshold = 0.75)
  
  names(clust) <- sapply(names(clust), function(x) nams[which(nams[,1]==x), 2])
  grp <- split(names(clust), clust)
  
  simu <- readRDS(sim)
  clus <- simu$gene_list
  seqs <- unlist(clus)
  nseqs <- length(seqs)
  mat_expect <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
  for (i in 1:length(clus)){
    ge <- clus[[i]]
    mat_expect[ge, ge] <- 1L
  }
  
  mat_observed <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
  for (i in 1:length(grp)){
    ge <- grp[[i]]
    mat_observed[ge, ge] <- 1L
  }
  
  
  res <- mat_expect %contingency% mat_observed
  return(res)
}


bench_micropan_blast_average <- function(din, dout = 'micropanBlast_vs_sim'){
  sim <- list.files(path = din, pattern = '[.]RDS$', full.names = TRUE)
  
  dd <- paste0(dout, '/', din)
  prep_faas <- list.files(dd, pattern = '[.]faa$', full.names = TRUE)
  nams <- lapply(prep_faas, function(x){
    fa <- read.fasta(x)
    an <- sapply(fa, attr, 'Annot')
    spl <- strsplit(sub('^>', '', an), ' ')
    mm <- do.call(rbind, spl)
    rownames(mm) <- NULL
    mm
  })
  
  nams <- do.call(rbind, nams)
  
  blout_dir <- paste(dd, 'blast_out', sep = '/')
  bl <- list.files(path = blout_dir, full.names = TRUE)
  df <- bDist(bl, verbose = FALSE)
  clust <- bClust(df, linkage = 'average', threshold = 0.75)
  
  names(clust) <- sapply(names(clust), function(x) nams[which(nams[,1]==x), 2])
  grp <- split(names(clust), clust)
  
  simu <- readRDS(sim)
  clus <- simu$gene_list
  seqs <- unlist(clus)
  nseqs <- length(seqs)
  mat_expect <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
  for (i in 1:length(clus)){
    ge <- clus[[i]]
    mat_expect[ge, ge] <- 1L
  }
  
  mat_observed <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
  for (i in 1:length(grp)){
    ge <- grp[[i]]
    mat_observed[ge, ge] <- 1L
  }
  
  
  res <- mat_expect %contingency% mat_observed
  return(res)
}



bench_micropan_hmmer <- function(din, dout = 'micropanHmmer_vs_sim', hmm_pfam){
  
  require(micropan)
  require(seqinr)
  
  faas <- list.files(path = din, pattern = '[.]faa$', full.names = TRUE)
  sim <- list.files(path = din, pattern = '[.]RDS$', full.names = TRUE)
  
  if (!dir.exists(dout)){
    dir.create(dout)
  }
  # dn <- dirname(gffs[1])
  dd <- paste0(dout, '/', din)
  if (dir.exists(dd)) unlink(dd, recursive = TRUE)
  dir.create(dd)
  
  
  gid <- paste0("%0",nchar(length(faas)),"d")
  gid <- paste0('GID',sprintf(gid,1:length(faas)))
  ref <- cbind(basename(faas), gid)
  write.csv(ref, file = paste0(dd, '/ref_gid.csv'),quote = F)
  
  prep_faas <- sapply(faas, function(x){
    ff <- basename(x)
    gd <- ref[which(ref[, 1] == ff), 2]
    out <- paste(dd, ff, sep = '/')
    panPrep(in.file = x, GID.tag = gd, out.file = out)
  })
  
  nams <- lapply(prep_faas, function(x){
    fa <- read.fasta(x)
    an <- sapply(fa, attr, 'Annot')
    spl <- strsplit(sub('^>', '', an), ' ')
    mm <- do.call(rbind, spl)
    rownames(mm) <- NULL
    mm
  })
  
  nams <- do.call(rbind, nams)
  
  hmmout_dir <- paste(dd, 'hmmer_out', sep = '/')
  if (dir.exists(hmmout_dir)) unlink(hmmout_dir, recursive = TRUE)
  dir.create(hmmout_dir)
  
  runMicropanHmmer <- function(prep_faas, oblf, db){
    hmmerScan(in.files = prep_faas, 
              db = db, 
              out.folder = oblf, 
              threads = 4L,
              verbose = FALSE)
    pfam.files <- list.files(path = oblf, full.names = TRUE)
    pfam.table <- lapply(pfam.files, function(i){
      tab <- readHmmer(i)
      tab <- hmmerCleanOverlap(tab)
      tab
    })
    pfam.table <- do.call(rbind, pfam.table)
    
    cluster.pfam <- dClust(pfam.table)
    return(cluster.pfam)
  }
  
  clust <- runMicropanHmmer(prep_faas = prep_faas, oblf = hmmout_dir, db = hmm_pfam)
  names(clust) <- sapply(names(clust), function(x) nams[which(nams[,1]==x), 2])
  grp <- split(names(clust), clust)
  
  simu <- readRDS(sim)
  clus <- simu$gene_list
  seqs <- unlist(clus)
  nseqs <- length(seqs)
  mat_expect <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
  for (i in 1:length(clus)){
    ge <- clus[[i]]
    mat_expect[ge, ge] <- 1L
  }
  
  mat_observed <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
  for (i in 1:length(grp)){
    ge <- grp[[i]]
    mat_observed[ge, ge] <- 1L
  }
  
  
  res <- mat_expect %contingency% mat_observed
  return(res)
  
}


# bench_FindMyFriends_KmerSplit <- function(din, dout = 'FindMyFriends_vs_sim'){
#   
#   require(FindMyFriends)
#   require(seqinr)
#   
#   faas <- list.files(path = din, pattern = '[.]faa$', full.names = TRUE)
#   gffs <- list.files(path = din, pattern = '[.]gff$', full.names = TRUE)
#   sim <- list.files(path = din, pattern = '[.]RDS$', full.names = TRUE)
#   
#   if (!dir.exists(dout)){
#     dir.create(dout)
#   }
#   # dn <- dirname(gffs[1])
#   dd <- paste0(dout, '/', din)
#   if (dir.exists(dd)) unlink(dd, recursive = TRUE)
#   dir.create(dd)
#   
#   file.copy(faas, dd, overwrite = TRUE)
#   faas2 <- list.files(path = dd, pattern = '[.]faa$', full.names = TRUE)
#   # Format headers
#   nams <- lapply(faas2, function(x){
#     gf <- grep(sub('[.]faa', '.gff', basename(x)), gffs, value = TRUE)
#     rl <- readLines(gf)
#     df <- pewit:::extractGffTable(rl)
#     nhe <- paste0(df$ID,'__',df$Contig, ' # ', df$From, ' # ', df$To, ' # ', ifelse(df$Strand=='+', 1, -1))
#     ff <- read.fasta(x)
#     write.fasta(ff, names = nhe, file.out = x)
#     cbind(nhe, df$ID)
#   })
#   nams <- do.call(rbind, nams)
#   
#   run_fmf <- function(faas){
#     pan <- FindMyFriends::pangenome(faas, 
#                                     translated = TRUE, 
#                                     geneLocation = 'prodigal', 
#                                     lowMem = F)
#     cdh <- cdhitGrouping(pan, kmerSize = 5)
#     nspl <- kmerSplit(cdh)
#     nspl
#   }
#   
#   clust <- run_fmf(faas = faas2)
#   grp <- lapply(genes(clust, split='group'), names)
#   grp <- lapply(grp, function(x) sub('__.+\\d', '', x))
#   
#   simu <- readRDS(sim)
#   clus <- simu$gene_list
#   seqs <- unlist(clus)
#   nseqs <- length(seqs)
#   mat_expect <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
#   for (i in 1:length(clus)){
#     ge <- clus[[i]]
#     mat_expect[ge, ge] <- 1L
#   }
#   
#   mat_observed <- matrix(0L, nrow = nseqs, ncol = nseqs, dimnames = list(seqs, seqs))
#   for (i in 1:length(grp)){
#     ge <- grp[[i]]
#     mat_observed[ge, ge] <- 1L
#   }
#   
#   
#   res <- mat_expect %contingency% mat_observed
#   return(res)
#   
# }


