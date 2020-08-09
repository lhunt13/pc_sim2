source('functions_gcomp.R')
library(data.table)
library(speedglm)
library(getopt)
nsim = 20000
B = 1000
n = 3000
K = 90
truth_L = readRDS('truth_L.rds')
truth_A = readRDS('truth_A.rds')

# set seed
args <- commandArgs(trailingOnly = TRUE)
boot.index <- as.numeric(args[1])
set.seed(boot.index)

# simulate data
data_L = sim_data(n,truth_L,K)
data_A = sim_data(n,truth_A,K)

# perform analysis
results = do_analysis(data_L,data_A,nsim,K)

# perform bootstrap
system.time({
  boot_results = matrix(NA,nrow=B,ncol=12)
  for(b in 1:B){
    # simulate data
    set.seed(10000 + b)
    boot_data_L = sim_data(n,results$fits_L,K)
    boot_data_A = sim_data(n,results$fits_A,K)
    
    boot_results[b,] = do_analysis(boot_data_L,boot_data_A,nsim,K)$estimates
  }
})

# get confidence intervals
bounds = apply(boot_results,2,function(x){quantile(x,probs = c(.025,.975))})
out = rbind(bounds,results$estimates)
out = c(out)
names(out) = paste0(c('L','U',''),rep(names(results$estimates),each=3))

# store out
# save file in "out" directory with file name "run-<boot.index>.rds"
out.file <- file.path("out", paste0("run-", boot.index, ".rds"))
saveRDS(out, out.file)
# quit R and don't save workspace
quit('no')



