
#r = getOption("repos")
#r["CRAN"] = "http://cran.us.r-project.org"
#options(repos = r)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("LEA")

library(LEA)

### Example of analysis using lfmm2 ###
args <- commandArgs(trailingOnly = TRUE)
print(args)

#mat <- as.matrix(read.csv(file=args[1], header = FALSE, sep=" "))
mat_pruned <-as.matrix(read.table(file=args[1], header=FALSE, sep=" "))
env <- as.matrix(read.table(file=args[2], header = FALSE, sep=" "))

class(mat_pruned) <- "integer"
class(env) <- "numeric"

#print(dim(mat_pruned))
# Apply the function to each row of the matrix
#apply(mat_pruned, 1, print_first_and_last)


######################################
# Fitting an LFMM with K factors #
######################################

#we really only see 3 maybe 4 actual clusters in the PCA. Screeplot suggests 3.  
#pass K as a paramter in the future, 6 is maybe too high, using tracy-widom test etc. 
#if we let X be the response this is still fine, since we aren't actually fitting anything, just computing the latent factors.
mod2 <- lfmm2(input = mat_pruned, env = env, K = 6)

saveRDS(mod2, args[3])
