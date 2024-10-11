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

#suppress warnings for now
options(warn=-1)
#options(warn=1)

mod2 <-readRDS(args[1])
mat <-as.matrix(read.table(file=args[2], header=FALSE, sep=" "))
env <- as.matrix(read.table(file=args[3], header = FALSE, sep=" "))

class(mat) <- "numeric"
class(env) <- "numeric"

###########################################
# Compute p-values given K latent factors #
###########################################

#default family is Binomial(logit), so we normalized the genotype counts by 2 since y needs to be between 0-1 
#run on pruned for testing
#may want to normalize genotypes differently? 


#pv <- lfmm2.test(object = mod2, input = mat/2.0, env = env, linear = FALSE, genomic.control = TRUE, family=quasibinomial(link="logit") )
#LINEAR TEST lfmm2 preprint did a case-control, because of course our Y matrix is the genotypes.  
#pv <- lfmm2.test(object = mod2, input = mat, env = env, linear = TRUE, genomic.control = TRUE)

#switching env and mat so we can do normal logistic regression. Shouldn't change results much at all but we can get effect sizes.
#pv <- lfmm2.test(object = mod2, input = env, env = mat, linear = FALSE, genomic.control = TRUE, family=quasibinomial(link="logit"))

genomic_control <- TRUE

#switching env and mat for logistic regression 
p <- ncol(mat)
p_value <- NULL
z_score <- NULL
d <- ncol(env) #number of env variables which is just 1

for (j in 1:p) {
	mod_glm <- glm(env ~., data = data.frame(mat[, j], mod2@U), family=binomial(link="logit"))
	sm <- summary(mod_glm)
#	print(sm)
	p_value <- rbind(p_value, sm$coeff[2:(d + 1), 4])
	z_score <- rbind(z_score, sm$coeff[2:(d +1), 3])

	#print(p_value)
}

if (genomic_control) {
	#if (d == 1) {
	#always 1 for us 
	gif <- median(z_score^2)/qchisq(0.5, df = 1, lower.tail = FALSE)
	#}
	p_value <- pchisq(z_score^2/gif, df = 1, lower.tail = FALSE)

}
res <- list(pvalues = p_value, zscores = z_score) 
write.table(res$pvalues, file=args[4], row.names=FALSE, col.names=FALSE)
