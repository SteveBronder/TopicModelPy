library(ggplot2)
library(ggthemes)
library(xtable)
library(coda)
library(data.table)

citations <- read.csv("./citations.csv")
likelihood <- read.csv("./likelihood.csv")
dtm_dat <- read.csv("./dtm_dat.csv")
phi <- read.csv("./phi_new.csv")
phi_old <- read.csv("./phi.csv")
theta <- read.csv("./theta.csv")
doc_length <- apply(dtm_dat[,-1],1,sum)
vocab <- colnames(phi)[-1]
term_frequency <- apply(dtm_dat[,-1],2,sum)

# Fix rounding error
phi <- phi[,-1]/apply(phi[,-1],1,sum)
phi_old <- phi_old[,-1]/apply(phi_old[,-1],2,sum)

#make and merge citations
citations[,2]<-as.character(citations[,2])
colnames(citations)[2] <- "source"
dtm_dat[,1] <- as.character(dtm_dat[,1])
dtm_dat[,1] <- substr(dtm_dat[,1],35,49)
dtm_dat[,1] <- gsub("_","/",dtm_dat[,1])

dtm_dat2 <- data.table(dtm_dat)
citations2 <- data.table(citations)
setkey(dtm_dat2,"source")
setkey(citations2,"source")
dtm_cit<-citations2[,c("source","title"),with=FALSE][dtm_dat2]
theta$X <- dtm_cit[,"title",with=FALSE]

ggplot(likelihood, aes(x=X, y=X0)) + geom_line() + theme_tufte() + xlab("Iterations") + ylab("Likelihood")

theta_top <- list()
for (i in 2:6){
  theta_test <- theta[order(theta[,i],decreasing=TRUE),"X"]
  theta_test <- as.data.frame(head(theta_test[title != "",],n=10))
  theta_top[I(i-1)] <- data.frame("topic_1" = as.character(theta_test$title),stringsAsFactors=FALSE)
}
lapply(1:5,function(i) (xtable(data.frame(theta_top[[i]]))))
theta_top <-do.call(cbind,theta_top)

xtable(theta_top)
# Geweke score
geweke.diag(mcmc(likelihood$X0))
heidel.diag(mcmc(likelihood$X0))

t_phio <- t(phi_old)

phi_top <- list()
for (i in 1:5){
  phi_top[i] <- as.data.frame( rownames(head(t_phio[order(t_phio[,i],decreasing=TRUE),],n=10)),stringAsFactors = FALSE)
}

phi_top <- lapply(phi_top, as.character)
phi_top <- do.call(cbind,phi_top)
phi_top <- as.data.frame(phi_top)
colnames(phi_top) <- paste0("topic_",1:5)

xtable(phi_top)

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
