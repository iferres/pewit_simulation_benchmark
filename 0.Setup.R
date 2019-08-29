################
## Containers ##
################

container_dir <- 'singularity_containers'
dir.create(container_dir)

## Repositories ##

# gff3toembl
gff3toembl <- list(repo = 'docker://sangerpathogens/gff3toembl', 
                   name = paste(container_dir, 'gfftoembl.sif', sep = '/'))
# roary
roary <- list(repo = 'shub://iferres/Singularity_recipes:roary',
              name = paste(container_dir, 'roary.sif', sep = '/'))
# panx
panx <- list(repo = 'shub://iferres/Singularity_recipes:panx',
             name = paste(container_dir, 'panx.sif', sep = '/'))

simg <- list(gff3toembl, roary, panx)

## Pull containers ##

lapply(simg, function(x){
  run <- paste('singularity pull', x$name, x$repo)
  system(run)
})



################
## R packages ##
################

## TO TEST ##

# PEWIT
if (!require(devtools)){
  install.packages(devtools)
  library(devtools)
}
install_github('iferres/pewit', ref = 'bcda1ad')


# MICROPAN
install.packages('micropan')


# FindMyFriends (dev)
devtools::install_github('thomasp85/FindMyFriends')



## Other packages needed ##
install.packages('seqinr')