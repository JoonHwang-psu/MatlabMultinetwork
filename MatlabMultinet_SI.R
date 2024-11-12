##### R codes for supplementary information (241001) of Matlab Multilayered network study


##### SI Table 1. Descriptive statistics of household-level variables
descript_node <- node_info[,c("HH_size", "head_age", "WEALTH_TOTAL_VAL_USD", "land", "salaried", "status")]

stargazer(descript_node, summary=TRUE, digits=2, 
          covariate.labels=c("Household size", "Household head age", "Household wealth (USD)", "Land ownership", "Salaried income", "Social status"),
          out="BD_ERGM_summary.htm")

##### SI Table 2. Network characteristics of the households depending on socioeconomic status
finsupp <- read.csv("https://raw.githubusercontent.com/JoonHwang-psu/MatlabMultinetwork/refs/heads/main/BD_finsupnet_summary.csv", header=TRUE, fileEncoding="UTF-8-BOM")
matsupp <- read.csv("https://raw.githubusercontent.com/JoonHwang-psu/MatlabMultinetwork/refs/heads/main/BD_matsupnet_summary.csv", header=TRUE, fileEncoding="UTF-8-BOM")
labsupp <- read.csv("https://raw.githubusercontent.com/JoonHwang-psu/MatlabMultinetwork/refs/heads/main/BD_labsupnet_summary.csv", header=TRUE, fileEncoding="UTF-8-BOM")
### these summaries were created using Gephi

## network characteristics in financial support layer

### land ownership
landfin<-by(finsupp[, c("indegree", "outdegree", "eigencentrality", "betweenesscentrality", "Hub", "Authority")], finsupp$land, colMeans)
(landfin <- do.call(rbind, landfin))
wilcox.test(indegree~land, data=finsupp)
wilcox.test(outdegree~land, data=finsupp)
wilcox.test(eigencentrality~land, data=finsupp)
wilcox.test(betweenesscentrality~land, data=finsupp)
wilcox.test(Hub~land, data=finsupp)
wilcox.test(Authority~land, data=finsupp)

### salaried income
salfin<-by(finsupp[, c("indegree", "outdegree", "eigencentrality", "betweenesscentrality", "Hub", "Authority")], finsupp$salaried, colMeans)
(salfin <- do.call(rbind, salfin))
wilcox.test(indegree~salaried, data=finsupp)
wilcox.test(outdegree~salaried, data=finsupp)
wilcox.test(eigencentrality~salaried, data=finsupp)
wilcox.test(betweenesscentrality~salaried, data=finsupp)
wilcox.test(Hub~salaried, data=finsupp)
wilcox.test(Authority~salaried, data=finsupp)

### institutional social status
statfin<-by(finsupp[, c("indegree", "outdegree", "eigencentrality", "betweenesscentrality", "Hub", "Authority")], finsupp$status, colMeans)
(statfin <- do.call(rbind, statfin))
wilcox.test(indegree~status, data=finsupp)
wilcox.test(outdegree~status, data=finsupp)
wilcox.test(eigencentrality~status, data=finsupp)
wilcox.test(betweenesscentrality~status, data=finsupp)
wilcox.test(Hub~status, data=finsupp)
wilcox.test(Authority~status, data=finsupp)

## network characteristics in material support layer

### land ownership
landmat<-by(matsupp[, c("indegree", "outdegree", "eigencentrality", "betweenesscentrality", "Hub", "Authority")], matsupp$land, colMeans)
(landmat <- do.call(rbind, landmat))
wilcox.test(indegree~land, data=matsupp)
wilcox.test(outdegree~land, data=matsupp)
wilcox.test(eigencentrality~land, data=matsupp)
wilcox.test(betweenesscentrality~land, data=matsupp)
wilcox.test(Hub~land, data=matsupp)
wilcox.test(Authority~land, data=matsupp)

### salaried income
salmat<-by(matsupp[, c("indegree", "outdegree", "eigencentrality", "betweenesscentrality", "Hub", "Authority")], matsupp$salaried, colMeans)
(salmat <- do.call(rbind, salmat))
wilcox.test(indegree~salaried, data=matsupp)
wilcox.test(outdegree~salaried, data=matsupp)
wilcox.test(eigencentrality~salaried, data=matsupp)
wilcox.test(betweenesscentrality~salaried, data=matsupp)
wilcox.test(Hub~salaried, data=matsupp)
wilcox.test(Authority~salaried, data=matsupp)

### institutional social status
statmat<-by(matsupp[, c("indegree", "outdegree", "eigencentrality", "betweenesscentrality", "Hub", "Authority")], matsupp$status, colMeans)
(statmat <- do.call(rbind, statmat))
wilcox.test(indegree~status, data=matsupp)
wilcox.test(outdegree~status, data=matsupp)
wilcox.test(eigencentrality~status, data=matsupp)
wilcox.test(betweenesscentrality~status, data=matsupp)
wilcox.test(Hub~status, data=matsupp)
wilcox.test(Authority~status, data=matsupp)

## network characteristics in labor support layer

### land ownership
landlab<-by(labsupp[, c("indegree", "outdegree", "eigencentrality", "betweenesscentrality", "Hub", "Authority")], labsupp$land, colMeans)
(landlab <- do.call(rbind, landlab))
wilcox.test(indegree~land, data=labsupp)
wilcox.test(outdegree~land, data=labsupp)
wilcox.test(eigencentrality~land, data=labsupp)
wilcox.test(betweenesscentrality~land, data=labsupp)
wilcox.test(Hub~land, data=labsupp)
wilcox.test(Authority~land, data=labsupp)

### salaried income
sallab<-by(labsupp[, c("indegree", "outdegree", "eigencentrality", "betweenesscentrality", "Hub", "Authority")], labsupp$salaried, colMeans)
(sallab <- do.call(rbind, sallab))
wilcox.test(indegree~salaried, data=labsupp)
wilcox.test(outdegree~salaried, data=labsupp)
wilcox.test(eigencentrality~salaried, data=labsupp)
wilcox.test(betweenesscentrality~salaried, data=labsupp)
wilcox.test(Hub~salaried, data=labsupp)
wilcox.test(Authority~salaried, data=labsupp)

### institutional social status
statlab<-by(labsupp[, c("indegree", "outdegree", "eigencentrality", "betweenesscentrality", "Hub", "Authority")], labsupp$status, colMeans)
(statlab <- do.call(rbind, statlab))
wilcox.test(indegree~status, data=labsupp)
wilcox.test(outdegree~status, data=labsupp)
wilcox.test(eigencentrality~status, data=labsupp)
wilcox.test(betweenesscentrality~status, data=labsupp)
wilcox.test(Hub~status, data=labsupp)
wilcox.test(Authority~status, data=labsupp)


##### SI Table 3. connections to non-villagers and institutions

external <- read.csv("https://raw.githubusercontent.com/JoonHwang-psu/MatlabMultinetwork/refs/heads/main/BD_external_summary.csv", header=TRUE, fileEncoding="UTF-8-BOM")

landext <- by(external[,c("outside", "institution")], external$land, colMeans)
(landext <- do.call(rbind, landext))
wilcox.test(outside~land, data=external)
wilcox.test(institution~land, data=external)

salext <- by(external[,c("outside", "institution")], external$salaried, colMeans)
(salext <- do.call(rbind, salext))
wilcox.test(outside~salaried, data=external)
wilcox.test(institution~salaried, data=external)

statext <- by(external[,c("outside", "institution")], external$status, colMeans)
(statext <- do.call(rbind, statext))
wilcox.test(outside~status, data=external)
wilcox.test(institution~status, data=external)

##### SI Figure 3-6 & SI Table 5-8. goodness-of-fit test for Matlab multilayered support network analysis (240925)

##### To run the gof test, I recommend keeping Rdata from the main analysis,
##### since the gof uses the fundamental model results to check whether randomly simulated networks structurally resemble empirical network

set.seed(NULL)

nsim<-100

sim.gof.model <- simulate(fundamental, nsim=nsim, output="network")

# extracting single layer networks from simulated multilayered networks

sim.cash = list()
sim.material = list()
sim.help = list()

for (i in 1:nsim) {
  type <- sim.gof.model[[i]] %v% ".LayerName"
  sim.cash[[i]] <- sim.gof.model[[i]] %s% which(type=="cash")
  sim.material[[i]] <- sim.gof.model[[i]] %s% which(type=="material")
  sim.help[[i]] <- sim.gof.model[[i]] %s% which(type=="help")
}

sim.cash.net <- network.list(sim.cash)
sim.material.net <- network.list(sim.material)
sim.help.net <- network.list(sim.help)

# gof test for in- and out-degrees

sim.cash.max.indeg <- numeric()
sim.material.max.indeg <- numeric()
sim.help.max.indeg <- numeric()

sim.cash.max.odeg <- numeric()
sim.material.max.odeg <- numeric()
sim.help.max.odeg <- numeric()

for (i in 1:nsim) {
  sim.cash.max.indeg[i] <- max(degree(sim.cash[[i]], cmode="indegree"))
  sim.material.max.indeg[i] <- max(degree(sim.material[[i]], cmode="indegree"))
  sim.help.max.indeg[i] <- max(degree(sim.help[[i]], cmode="indegree"))
  sim.cash.max.odeg[i] <- max(degree(sim.cash[[i]], cmode="outdegree"))
  sim.material.max.odeg[i] <- max(degree(sim.material[[i]], cmode="outdegree"))
  sim.help.max.odeg[i] <- max(degree(sim.help[[i]], cmode="outdegree"))
}

emp.cash.max.indeg <- max(degree(cash_network, cmode="indegree"))
emp.material.max.indeg <- max(degree(material_network, cmode="indegree"))
emp.help.max.indeg <- max(degree(help_network, cmode="indegree"))

emp.cash.max.odeg <- max(degree(cash_network, cmode="outdegree"))
emp.material.max.odeg <- max(degree(material_network, cmode="outdegree"))
emp.help.max.odeg <- max(degree(help_network, cmode="outdegree"))

cash.max.indeg <- max(sim.cash.max.indeg, emp.cash.max.indeg)
material.max.indeg <- max(sim.material.max.indeg, emp.material.max.indeg)
help.max.indeg <- max(sim.help.max.indeg, emp.help.max.indeg)

cash.max.odeg <- max(sim.cash.max.odeg, emp.cash.max.odeg)
material.max.odeg <- max(sim.material.max.odeg, emp.material.max.odeg)
help.max.odeg <- max(sim.help.max.odeg, emp.help.max.odeg)

sim.cash.indeg <- summary(sim.cash.net ~ idegree(0:cash.max.indeg))
sim.material.indeg <- summary(sim.material.net ~ idegree(0:material.max.indeg))
sim.help.indeg <- summary(sim.help.net ~ idegree(0:help.max.indeg))

emp.cash.indeg <- summary(cash_network ~ idegree(0:cash.max.indeg))
emp.material.indeg <- summary(material_network ~ idegree(0:material.max.indeg))
emp.help.indeg <- summary(help_network ~ idegree(0:help.max.indeg))

sim.cash.odeg <- summary(sim.cash.net ~ odegree(0:cash.max.odeg))
sim.material.odeg <- summary(sim.material.net ~ odegree(0:material.max.odeg))
sim.help.odeg <- summary(sim.help.net ~ odegree(0:help.max.odeg))

emp.cash.odeg <- summary(cash_network ~ odegree(0:cash.max.odeg))
emp.material.odeg <- summary(material_network ~ odegree(0:material.max.odeg))
emp.help.odeg <- summary(help_network ~ odegree(0:help.max.odeg))

# empirical p-value testing for indegree

labelsize <- 1.5
par(mfrow=c(3,1), mar=c(5,3,1,1))

pval <- apply(sim.cash.indeg <= emp.cash.indeg[col(sim.cash.indeg)],2, mean)
pval.top <- apply(sim.cash.indeg >= emp.cash.indeg[col(sim.cash.indeg)],2, mean)
pval <- cbind(emp.cash.indeg,apply(sim.cash.indeg, 2,min), apply(sim.cash.indeg, 2,mean),
              apply(sim.cash.indeg, 2,max), pmin(1,2*pmin(pval,pval.top)))
dimnames(pval)[[2]] <- c("obs","min","mean","max","MC p-value")
pval

psim <- sweep(sim.cash.indeg,1,apply(sim.cash.indeg,1,sum),"/")
pobs <- emp.cash.indeg/sum(emp.cash.indeg)
bds <- apply(psim,2,quantile,probs=c(0.025,0.975))

boxplot(psim, outline=FALSE, xlab="Financial support network", cex.axis = labelsize, cex.lab = labelsize)
points(bds[1,], pch=1, cex=0.75)
points(bds[2,], pch=1, cex=0.75)
lines(bds[1,], pch = 18,lty=1,lwd=1,col="gray75")
lines(bds[2,], pch = 18,lty=1,lwd=1,col="gray75")
points(pobs, pch = 16,cex=0.75)
lines(pobs, lty = 1,lwd=3)
points(colMeans(psim), pch=18, cex=2, col="blue")

pval <- apply(sim.material.indeg <= emp.material.indeg[col(sim.material.indeg)],2, mean)
pval.top <- apply(sim.material.indeg >= emp.material.indeg[col(sim.material.indeg)],2, mean)
pval<- cbind(emp.material.indeg,apply(sim.material.indeg, 2,min), apply(sim.material.indeg, 2,mean),
                      apply(sim.material.indeg, 2,max), pmin(1,2*pmin(pval,pval.top)))
dimnames(pval)[[2]] <- c("obs","min","mean","max","MC p-value")
pval

psim <- sweep(sim.material.indeg,1,apply(sim.material.indeg,1,sum),"/")
pobs <- emp.material.indeg/sum(emp.material.indeg)
bds <- apply(psim,2,quantile,probs=c(0.025,0.975))

boxplot(psim, outline=FALSE, xlab="Material support network", cex.axis = labelsize, cex.lab = labelsize)
points(bds[1,], pch=1, cex=0.75)
points(bds[2,], pch=1, cex=0.75)
lines(bds[1,], pch = 18,lty=1,lwd=1,col="gray75")
lines(bds[2,], pch = 18,lty=1,lwd=1,col="gray75")
points(pobs, pch = 16,cex=0.75)
lines(pobs, lty = 1,lwd=3)
points(colMeans(psim), pch=18, cex=2, col="blue")

pval <- apply(sim.help.indeg <= emp.help.indeg[col(sim.help.indeg)],2, mean)
pval.top <- apply(sim.help.indeg >= emp.help.indeg[col(sim.help.indeg)],2, mean)
pval <- cbind(emp.help.indeg,apply(sim.help.indeg, 2,min), apply(sim.help.indeg, 2,mean),
                      apply(sim.help.indeg, 2,max), pmin(1,2*pmin(pval,pval.top)))
dimnames(pval)[[2]] <- c("obs","min","mean","max","MC p-value")
pval

psim <- sweep(sim.help.indeg,1,apply(sim.help.indeg,1,sum),"/")
pobs <- emp.help.indeg/sum(emp.help.indeg)
bds <- apply(psim,2,quantile,probs=c(0.025,0.975))

boxplot(psim, outline=FALSE, xlab="Labor support network", cex.axis = labelsize, cex.lab = labelsize)
points(bds[1,], pch=1, cex=0.75)
points(bds[2,], pch=1, cex=0.75)
lines(bds[1,], pch = 18,lty=1,lwd=1,col="gray75")
lines(bds[2,], pch = 18,lty=1,lwd=1,col="gray75")
points(pobs, pch = 16,cex=0.75)
lines(pobs, lty = 1,lwd=3)
points(colMeans(psim), pch=18, cex=2, col="blue")

# empirical p-value testing for outdegrees

pval <- apply(sim.cash.odeg <= emp.cash.odeg[col(sim.cash.odeg)],2, mean)
pval.top <- apply(sim.cash.odeg >= emp.cash.odeg[col(sim.cash.odeg)],2, mean)
pval <- cbind(emp.cash.odeg,apply(sim.cash.odeg, 2,min), apply(sim.cash.odeg, 2,mean),
                      apply(sim.cash.odeg, 2,max), pmin(1,2*pmin(pval,pval.top)))
dimnames(pval)[[2]] <- c("obs","min","mean","max","MC p-value")
pval

psim <- sweep(sim.cash.odeg,1,apply(sim.cash.odeg,1,sum),"/")
pobs <- emp.cash.odeg/sum(emp.cash.odeg)
bds <- apply(psim,2,quantile,probs=c(0.025,0.975))

boxplot(psim, outline=FALSE, xlab="Financial support network", cex.axis = labelsize, cex.lab = labelsize)
points(bds[1,], pch=1, cex=0.75)
points(bds[2,], pch=1, cex=0.75)
lines(bds[1,], pch = 18,lty=1,lwd=1,col="gray75")
lines(bds[2,], pch = 18,lty=1,lwd=1,col="gray75")
points(pobs, pch = 16,cex=0.75)
lines(pobs, lty = 1,lwd=3)
points(colMeans(psim), pch=18, cex=2, col="blue")

pval <- apply(sim.material.odeg <= emp.material.odeg[col(sim.material.odeg)],2, mean)
pval.top <- apply(sim.material.odeg >= emp.material.odeg[col(sim.material.odeg)],2, mean)
pval <- cbind(emp.material.odeg,apply(sim.material.odeg, 2,min), apply(sim.material.odeg, 2,mean),
                      apply(sim.material.odeg, 2,max), pmin(1,2*pmin(pval,pval.top)))
dimnames(pval)[[2]] <- c("obs","min","mean","max","MC p-value")
pval

psim <- sweep(sim.material.odeg,1,apply(sim.material.odeg,1,sum),"/")
pobs <- emp.material.odeg/sum(emp.material.odeg)
bds <- apply(psim,2,quantile,probs=c(0.025,0.975))

boxplot(psim, outline=FALSE, xlab="Material support network", cex.axis = labelsize, cex.lab = labelsize)
points(bds[1,], pch=1, cex=0.75)
points(bds[2,], pch=1, cex=0.75)
lines(bds[1,], pch = 18,lty=1,lwd=1,col="gray75")
lines(bds[2,], pch = 18,lty=1,lwd=1,col="gray75")
points(pobs, pch = 16,cex=0.75)
lines(pobs, lty = 1,lwd=3)
points(colMeans(psim), pch=18, cex=2, col="blue")

pval <- apply(sim.help.odeg <= emp.help.odeg[col(sim.help.odeg)],2, mean)
pval.top <- apply(sim.help.odeg >= emp.help.odeg[col(sim.help.odeg)],2, mean)
pval <- cbind(emp.help.odeg,apply(sim.help.odeg, 2,min), apply(sim.help.odeg, 2,mean),
                      apply(sim.help.odeg, 2,max), pmin(1,2*pmin(pval,pval.top)))
dimnames(pval)[[2]] <- c("obs","min","mean","max","MC p-value")
pval

psim <- sweep(sim.help.odeg,1,apply(sim.help.odeg,1,sum),"/")
pobs <- emp.help.odeg/sum(emp.help.odeg)
bds <- apply(psim,2,quantile,probs=c(0.025,0.975))

boxplot(psim, outline=FALSE, xlab="Labor support network", cex.axis = labelsize, cex.lab = labelsize)
points(bds[1,], pch=1, cex=0.75)
points(bds[2,], pch=1, cex=0.75)
lines(bds[1,], pch = 18,lty=1,lwd=1,col="gray75")
lines(bds[2,], pch = 18,lty=1,lwd=1,col="gray75")
points(pobs, pch = 16,cex=0.75)
lines(pobs, lty = 1,lwd=3)
points(colMeans(psim), pch=18, cex=2, col="blue")


##### gof test for geodesic distance

emp.cash.geodist <- ergm.geodistdist(cash_network)
emp.material.geodist <- ergm.geodistdist(material_network)
emp.help.geodist <- ergm.geodistdist(help_network)

(emp.cash.geodist.max <-length(emp.cash.geodist[emp.cash.geodist != 0])-1)
(emp.material.geodist.max <- length(emp.material.geodist[emp.material.geodist != 0])-1)
(emp.help.geodist.max<-length(emp.help.geodist[emp.help.geodist != 0])-1)

sim.cash.geodist <- list()
sim.material.geodist <- list()
sim.help.geodist <- list()

for (i in 1:nsim){
  sim.cash.geodist[[i]] <- ergm.geodistdist(sim.cash.net[[i]])  
  sim.material.geodist[[i]] <- ergm.geodistdist(sim.material.net[[i]])  
  sim.help.geodist[[i]] <- ergm.geodistdist(sim.help.net[[i]])  
}

sim.cash.geodist <- data.frame(sim.cash.geodist)
sim.material.geodist <- data.frame(sim.material.geodist)
sim.help.geodist <- data.frame(sim.help.geodist)

colnames(sim.cash.geodist)<-1:nsim
colnames(sim.material.geodist)<-1:nsim
colnames(sim.help.geodist)<-1:nsim

sim.cash.geodist <- t(sim.cash.geodist)
sim.material.geodist <- t(sim.material.geodist)
sim.help.geodist <- t(sim.help.geodist)

#sim.cash.geodist.max <- as.numeric(nrow(sim.cash.geodist[rowSums(sim.cash.geodist[])>0,]))-1
#sim.material.geodist.max <- as.numeric(nrow(sim.material.geodist[rowSums(sim.material.geodist[])>0,]))-1
#sim.help.geodist.max <- as.numeric(nrow(sim.help.geodist[rowSums(sim.help.geodist[])>0,]))-1

cash.geodist.max <- max(emp.cash.geodist.max #, sim.cash.geodist.max
                        )
material.geodist.max <- max(emp.material.geodist.max #, sim.material.geodist.max
                            )
help.geodist.max <- max(emp.help.geodist.max #, sim.help.geodist.max
                        )

sim.cash.geodist <- sim.cash.geodist[,c(1:cash.geodist.max, "Inf")]
sim.material.geodist <- sim.material.geodist[,c(1:material.geodist.max, "Inf")]
sim.help.geodist <- sim.help.geodist[,c(1:help.geodist.max, "Inf")]

emp.cash.geodist <- emp.cash.geodist[c(1:cash.geodist.max, "Inf")]
emp.material.geodist <- emp.material.geodist[c(1:material.geodist.max, "Inf")]
emp.help.geodist <- emp.help.geodist[c(1:help.geodist.max, "Inf")]

# empirical p-value testing for geodesic distance

pval <- apply(sim.cash.geodist <= emp.cash.geodist[col(sim.cash.geodist)],2, mean)
pval.top <- apply(sim.cash.geodist >= emp.cash.geodist[col(sim.cash.geodist)],2, mean)
pval <- cbind(emp.cash.geodist,apply(sim.cash.geodist, 2,min), apply(sim.cash.geodist, 2,mean),
                      apply(sim.cash.geodist, 2,max), pmin(1,2*pmin(pval,pval.top)))
dimnames(pval)[[2]] <- c("obs","min","mean","max","MC p-value")
pval

psim <- sweep(sim.cash.geodist,1,apply(sim.cash.geodist,1,sum),"/")
pobs <- emp.cash.geodist/sum(emp.cash.geodist)
bds <- apply(psim,2,quantile,probs=c(0.025,0.975))

boxplot(psim[,1:cash.geodist.max], outline=FALSE, xlab="Financial support network", cex.axis = labelsize, cex.lab = labelsize)
points(bds[1,1:cash.geodist.max], pch=1, cex=0.75)
points(bds[2,1:cash.geodist.max], pch=1, cex=0.75)
lines(bds[1,1:cash.geodist.max], pch = 18,lty=1,lwd=1,col="gray75")
lines(bds[2,1:cash.geodist.max], pch = 18,lty=1,lwd=1,col="gray75")
points(pobs[1:cash.geodist.max], pch = 16,cex=0.75)
lines(pobs[1:cash.geodist.max], lty = 1,lwd=3)
points(colMeans(psim[,1:cash.geodist.max]), pch=18, cex=2, col="blue")

pval <- apply(sim.material.geodist <= emp.material.geodist[col(sim.material.geodist)],2, mean)
pval.top <- apply(sim.material.geodist >= emp.material.geodist[col(sim.material.geodist)],2, mean)
pval <- cbind(emp.material.geodist,apply(sim.material.geodist, 2,min), apply(sim.material.geodist, 2,mean),
                      apply(sim.material.geodist, 2,max), pmin(1,2*pmin(pval,pval.top)))
dimnames(pval)[[2]] <- c("obs","min","mean","max","MC p-value")
pval

psim <- sweep(sim.material.geodist,1,apply(sim.material.geodist,1,sum),"/")
pobs <- emp.material.geodist/sum(emp.material.geodist)
bds <- apply(psim,2,quantile,probs=c(0.025,0.975))

boxplot(psim[,1:material.geodist.max], outline=FALSE, xlab="Material support network", cex.axis = labelsize, cex.lab = labelsize)
points(bds[1,1:material.geodist.max], pch=1, cex=0.75)
points(bds[2,1:material.geodist.max], pch=1, cex=0.75)
lines(bds[1,1:material.geodist.max], pch = 18,lty=1,lwd=1,col="gray75")
lines(bds[2,1:material.geodist.max], pch = 18,lty=1,lwd=1,col="gray75")
points(pobs[1:material.geodist.max], pch = 16,cex=0.75)
lines(pobs[1:material.geodist.max], lty = 1,lwd=3)
points(colMeans(psim[,1:material.geodist.max]), pch=18, cex=2, col="blue")

pval <- apply(sim.help.geodist <= emp.help.geodist[col(sim.help.geodist)],2, mean)
pval.top <- apply(sim.help.geodist >= emp.help.geodist[col(sim.help.geodist)],2, mean)
pval <- cbind(emp.help.geodist,apply(sim.help.geodist, 2,min), apply(sim.help.geodist, 2,mean),
                      apply(sim.help.geodist, 2,max), pmin(1,2*pmin(pval,pval.top)))
dimnames(pval)[[2]] <- c("obs","min","mean","max","MC p-value")
pval

psim <- sweep(sim.help.geodist,1,apply(sim.help.geodist,1,sum),"/")
pobs <- emp.help.geodist/sum(emp.help.geodist)
bds <- apply(psim,2,quantile,probs=c(0.025,0.975))

boxplot(psim[,1:help.geodist.max], outline=FALSE, xlab="Labor support network", cex.axis = labelsize, cex.lab = labelsize)
points(bds[1,1:help.geodist.max], pch=1, cex=0.75)
points(bds[2,1:help.geodist.max], pch=1, cex=0.75)
lines(bds[1,1:help.geodist.max], pch = 18,lty=1,lwd=1,col="gray75")
lines(bds[2,1:help.geodist.max], pch = 18,lty=1,lwd=1,col="gray75")
points(pobs[1:help.geodist.max], pch = 16,cex=0.75)
lines(pobs[1:help.geodist.max], lty = 1,lwd=3)
points(colMeans(psim[,1:help.geodist.max]), pch=18, cex=2, col="blue")

##### gof test for edgewise shared partners

emp.cash.esp <- summary(cash_network~esp(0:100))
emp.material.esp <- summary(material_network~esp(0:100))
emp.help.esp <- summary(help_network~esp(0:100))

(emp.cash.esp.max <-length(emp.cash.esp[emp.cash.esp != 0]))
(emp.material.esp.max <- length(emp.material.esp[emp.material.esp != 0]))
(emp.help.esp.max<-length(emp.help.esp[emp.help.esp != 0]))

sim.cash.esp <- summary(sim.cash.net~esp(0:100))
sim.material.esp <- summary(sim.material.net~esp(0:100))
sim.help.esp <- summary(sim.help.net~esp(0:100))

(sim.cash.esp.max <- as.numeric(ncol(sim.cash.esp[,colSums(sim.cash.esp[])>0])))
(sim.material.esp.max <- as.numeric(ncol(sim.material.esp[,colSums(sim.material.esp[])>0])))
(sim.help.esp.max <- as.numeric(ncol(sim.help.esp[,colSums(sim.help.esp[])>0])))

(cash.esp.max <- max(emp.cash.esp.max, sim.cash.esp.max))
(material.esp.max <- max(emp.material.esp.max, sim.material.esp.max))
(help.esp.max <- max(emp.help.esp.max, sim.help.esp.max))

sim.cash.esp <- sim.cash.esp[,c(1:cash.esp.max)]
sim.material.esp <- sim.material.esp[,c(1:material.esp.max)]
sim.help.esp <- sim.help.esp[,c(1:help.esp.max)]

emp.cash.esp <- emp.cash.esp[c(1:cash.esp.max)]
emp.material.esp <- emp.material.esp[c(1:material.esp.max)]
emp.help.esp <- emp.help.esp[c(1:help.esp.max)]

# empirical p-value testing for edgewise shared partners

pval <- apply(sim.cash.esp <= emp.cash.esp[col(sim.cash.esp)],2, mean)
pval.top <- apply(sim.cash.esp >= emp.cash.esp[col(sim.cash.esp)],2, mean)
pval <- cbind(emp.cash.esp,apply(sim.cash.esp, 2,min), apply(sim.cash.esp, 2,mean),
                      apply(sim.cash.esp, 2,max), pmin(1,2*pmin(pval,pval.top)))
dimnames(pval)[[2]] <- c("obs","min","mean","max","MC p-value")
pval

psim <- sweep(sim.cash.esp,1,apply(sim.cash.esp,1,sum),"/")
pobs <- emp.cash.esp/sum(emp.cash.esp)
bds <- apply(psim,2,quantile,probs=c(0.025,0.975))

boxplot(psim, outline=FALSE, xlab="Financial support network", cex.axis = labelsize, cex.lab = labelsize)
points(bds[1,], pch=1, cex=0.75)
points(bds[2,], pch=1, cex=0.75)
lines(bds[1,], pch = 18,lty=1,lwd=1,col="gray75")
lines(bds[2,], pch = 18,lty=1,lwd=1,col="gray75")
points(pobs, pch = 16,cex=0.75)
lines(pobs, lty = 1,lwd=3)
points(colMeans(psim), pch=18, cex=2, col="blue")

pval <- apply(sim.material.esp <= emp.material.esp[col(sim.material.esp)],2, mean)
pval.top <- apply(sim.material.esp >= emp.material.esp[col(sim.material.esp)],2, mean)
pval <- cbind(emp.material.esp,apply(sim.material.esp, 2,min), apply(sim.material.esp, 2,mean),
                      apply(sim.material.esp, 2,max), pmin(1,2*pmin(pval,pval.top)))
dimnames(pval)[[2]] <- c("obs","min","mean","max","MC p-value")
pval

psim <- sweep(sim.material.esp,1,apply(sim.material.esp,1,sum),"/")
pobs <- emp.material.esp/sum(emp.material.esp)
bds <- apply(psim,2,quantile,probs=c(0.025,0.975))

boxplot(psim, outline=FALSE, xlab="Material support network", cex.axis = labelsize, cex.lab = labelsize)
points(bds[1,], pch=1, cex=0.75)
points(bds[2,], pch=1, cex=0.75)
lines(bds[1,], pch = 18,lty=1,lwd=1,col="gray75")
lines(bds[2,], pch = 18,lty=1,lwd=1,col="gray75")
points(pobs, pch = 16,cex=0.75)
lines(pobs, lty = 1,lwd=3)
points(colMeans(psim), pch=18, cex=2, col="blue")

pval <- apply(sim.help.esp <= emp.help.esp[col(sim.help.esp)],2, mean)
pval.top <- apply(sim.help.esp >= emp.help.esp[col(sim.help.esp)],2, mean)
pval <- cbind(emp.help.esp,apply(sim.help.esp, 2,min), apply(sim.help.esp, 2,mean),
                      apply(sim.help.esp, 2,max), pmin(1,2*pmin(pval,pval.top)))
dimnames(pval)[[2]] <- c("obs","min","mean","max","MC p-value")
pval

psim <- sweep(sim.help.esp,1,apply(sim.help.esp,1,sum),"/")
pobs <- emp.help.esp/sum(emp.help.esp)
bds <- apply(psim,2,quantile,probs=c(0.025,0.975))

boxplot(psim, outline=FALSE, xlab="Labor support network", cex.axis = labelsize, cex.lab = labelsize)
points(bds[1,], pch=1, cex=0.75)
points(bds[2,], pch=1, cex=0.75)
lines(bds[1,], pch = 18,lty=1,lwd=1,col="gray75")
lines(bds[2,], pch = 18,lty=1,lwd=1,col="gray75")
points(pobs, pch = 16,cex=0.75)
lines(pobs, lty = 1,lwd=3)
points(colMeans(psim), pch=18, cex=2, col="blue")

# density test

density <- gden(sim.gof.model)
summary(density)
gden(multi_network)

density.cash <- gden(sim.cash.net)
summary(density.cash)
emp.density.cash <- gden(cash_network)

density.material <- gden(sim.material.net)
summary(density.material)
emp.density.material <- gden(material_network)

density.help <- gden(sim.help.net)
summary(density.help)
emp.density.help <- gden(help_network)

sim.density.df <- data.frame(density.cash, density.material, density.help)
emp.density <- c(emp.density.cash, emp.density.material, emp.density.help)

par(mfrow = c(1, 1), mar=c(5,1,1,1))

boxplot(sim.density.df, names=c("financial", "material", "labor"), cex.axis = labelsize, cex.lab = labelsize)
points(emp.density, pch=18, cex=2, col=2)
legend("topright", bty = "n", pch = 18, pt.cex = 1.5, pt.lwd = 2, col = 2, cex = 0.75*labelsize,
       legend = "Empirical")

