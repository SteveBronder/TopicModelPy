library(ggplot2)
library(ggthemes)
library(coda)

likelihood <- read.csv("./likelihood.csv")
dtm_dat <- read.csv("./dtm_dat.csv")
phi <- read.csv("./phi_new.csv")
phi_old <- read.csv("./phi.csv")
# Fix rounding error
phi <- phi[,-1]/apply(phi[,-1],1,sum)
phi_old <- phi_old[,-1]/apply(phi_old[,-1],2,sum)

theta <- read.csv("./theta.csv")
doc_length <- apply(dtm_dat[,-1],1,sum)
vocab <- colnames(phi)[-1]
term_frequency <- apply(dtm_dat[,-1],2,sum)
ggplot(likelihood, aes(x=X, y=X0)) + geom_line() + theme_tufte() + xlab("Iterations") + ylab("Likelihood")

# Geweke score
geweke.diag(mcmc(likelihood$X0))

t_phio <- t(phi_old)

phi_top <- list()
for (i in 1:5){
  phi_top[i] <- as.data.frame( rownames(head(t_phio[order(t_phio[,i],decreasing=TRUE),],n=10)),stringAsFactors = FALSE)
}

phi_top <- lapply(phi_top, as.character)
phi_top <- do.call(cbind,phi_top)
phi_top <- as.data.frame(phi_top)
colnames(phi_top) <- paste0("topic_",1:5)

# Fancy Graph
##########################################
jsPCA2 <- function(phi) {
  # first, we compute a pairwise distance between topic distributions
  # using a symmetric version of KL-divergence
  # http://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence
  jensenShannon <- function(x, y) {
    m <- 0.5*(x + y)
    0.5*sum(x*log(x/m)) + 0.5*sum(y*log(y/m))
  }
  dist.mat <- stats::dist(x = phi_js, method = "canberra")
  # then, we reduce the K by K proximity matrix down to K by 2 using PCA
  phi_js <<-phi
  pca.fit <- stats::cmdscale(dist.mat, k = 2)
  data.frame(x = pca.fit[,1], y = pca.fit[,2])
}

dat_test <- createJSON(theta = theta[-1], phi = phi, doc.length = doc_length,
                       vocab = vocab, term.frequency = term_frequency,mds.method=jsPCA2)

serVis(dat_test)

test_phi <- t(phi)
test_phi <- test_phi[order(test_phi[,1],decreasing = TRUE),]
