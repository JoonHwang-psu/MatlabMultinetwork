##### Matlab multilayered support network analysis (240925)

library(igraph)
library(sna)
library(statnet)
library(intergraph)
library(stargazer)
library(latticeExtra)
library(ergm.multi)
library(network)

# 1. setting up the data

## 1.1. loading node attribute data
node_info <- read.csv("https://raw.githubusercontent.com/JoonHwang-psu/MatlabMultinetwork/refs/heads/main/BD_SharingUnit_NHB.csv", header=TRUE, fileEncoding="UTF-8-BOM")
node_info$wealth1000 <- (node_info$WEALTH_TOTAL_VAL_USD)/1000

colnames(node_info)[1] <- "name"

## 1.2. creating kinship networks
kin_edgelist<-read.csv("https://raw.githubusercontent.com/JoonHwang-psu/MatlabMultinetwork/refs/heads/main/BD_Edgelist_kinship_NHB.csv", header=TRUE)
kin_graph <- graph.data.frame(kin_edgelist, directed=FALSE)
kin_graph <- kin_graph - vertex("BDSU019")
kin_graph <- simplify(kin_graph, remove.multiple=TRUE, remove.loops=TRUE,edge.attr.comb = igraph_opt("edge.attr.comb"))
is_simple(kin_graph)

kin_network <- asNetwork(kin_graph)

## 1.3. creating distance networks
distance_edgelist<-read.csv("https://raw.githubusercontent.com/JoonHwang-psu/MatlabMultinetwork/refs/heads/main/BD_Distance_NHB.csv", header=TRUE)
distance_graph <- graph.data.frame(distance_edgelist, directed=FALSE)
distance_graph <- distance_graph - vertex("BDSU019")
distance_graph <- simplify(distance_graph, remove.multiple=TRUE, remove.loops=TRUE,edge.attr.comb = "mean")

distance_network <- asNetwork(distance_graph)

## 1.4. creating multilayered support network

### financial support network
cash_edgelist<-read.csv("https://raw.githubusercontent.com/JoonHwang-psu/MatlabMultinetwork/refs/heads/main/BD_Edgelist_Q1Q2_NHB.csv", header=TRUE)
cash_graph <- graph.data.frame(cash_edgelist, directed=TRUE)
cash_graph <- cash_graph-vertex("BDSU019")+vertex("BDSU044")
cash_graph <- simplify(cash_graph, remove.multiple=TRUE, remove.loops=TRUE,edge.attr.comb = igraph_opt("edge.attr.comb"))
is_simple(cash_graph)

cash_edgelist2 <- as.data.frame(get.edgelist(cash_graph))
cash_edgelist2$cash<-1
cash_edgelist2$material<-0
cash_edgelist2$help<-0

### material support network
material_edgelist<-read.csv("https://raw.githubusercontent.com/JoonHwang-psu/MatlabMultinetwork/refs/heads/main/BD_Edgelist_Q3Q4_NHB.csv", header=TRUE)
material_graph <- graph.data.frame(material_edgelist, directed=TRUE)
material_graph <- material_graph - vertex("BDSU019")+vertex("BDSU044")
material_graph <- simplify(material_graph, remove.multiple=TRUE, remove.loops=TRUE,edge.attr.comb = igraph_opt("edge.attr.comb"))
is_simple(material_graph)

material_edgelist2<-as.data.frame(get.edgelist(material_graph))
material_edgelist2$cash<-0
material_edgelist2$material<-1
material_edgelist2$help<-0

### labor support network
inv_help_edgelist<-read.csv("https://raw.githubusercontent.com/JoonHwang-psu/MatlabMultinetwork/refs/heads/main/BD_Edgelist_Q5Q6_NHB.csv", header=TRUE) ## tie direction in original edgelist is in reverse to the resource flow
inv_help_graph <- graph.data.frame(inv_help_edgelist, directed=TRUE)
inv_help_graph <- inv_help_graph - vertex("BDSU019") 
help_graph <- reverse_edges(inv_help_graph) ## so here I reversed it so it can match actual flow of labor service 
help_graph <- simplify(help_graph, remove.multiple=TRUE, remove.loops=TRUE,edge.attr.comb = igraph_opt("edge.attr.comb"))
is_simple(help_graph)

help_edgelist2<-as.data.frame(get.edgelist(help_graph))
help_edgelist2$cash<-0
help_edgelist2$material<-0
help_edgelist2$help<-1

total_edgelist<-rbind(cash_edgelist2,material_edgelist2,help_edgelist2) ### creating combined edgelist including all three types of support relationships
total_graph<-graph.data.frame(total_edgelist, directed=TRUE) ### and converting it into combined graph object

### attaching nodal attributes to combined graph
V(total_graph)$name
V(total_graph)$wealth1000 <- sapply(V(total_graph)$name, function(x) node_info$wealth1000[node_info$name == x])
V(total_graph)$land <- sapply(V(total_graph)$name, function(x) node_info$land[node_info$name == x])
V(total_graph)$salaried <- sapply(V(total_graph)$name, function(x) node_info$salaried[node_info$name == x])
V(total_graph)$status<- sapply(V(total_graph)$name, function(x) node_info$status[node_info$name == x])

### now converting combined graph into combined networks which do not have layer elements yet.
total_network<-asNetwork(total_graph)

(cash_network<-network_view(total_network, "cash"))
(material_network<-network_view(total_network, "material"))
(help_network<-network_view(total_network, "help"))

### converting combined networks into multilayered networks by specifying the layers
multi_network<-Layer(total_network, c("cash", "material", "help"))

summary(multi_network)

# 3. multilayered ERGM and motif analysis

## 3.1. fitting fundamental model using multilayered ERGM
fundamental <- ergm(multi_network~
                    L(~edgecov(distance_network, "distance")+edgecov(kin_network), c(~cash, ~material, ~help))
                    +L(~edges+isolates+esp(1)+mutual+diff("wealth1000", pow=1, dir="t-h", sign.action="identity"), ~cash)
                    +L(~edges+isolates+esp(1)+mutual+diff("wealth1000", pow=1, dir="t-h", sign.action="identity"), ~material)
                    +L(~edges+isolates+esp(1)+mutual+diff("wealth1000", pow=1, dir="t-h", sign.action="identity"), ~help)
                    +L(~edges, ~material&cash)
                    +L(~edges, ~help&cash)
                    +L(~edges, ~help&material)
                    +mutualL(Ls=c(~material, ~cash))
                    +mutualL(Ls=c(~help, ~cash))
                    +mutualL(Ls=c(~help, ~material)))
summary(fundamental)

round(exp(cbind(coef(fundamental), confint(fundamental))), digits=3)

## 3.2. simulating random networks based on the fundamental model
set.seed(890716) # to see slightly different results, try set.seed(NULL)

sim.model<-simulate(fundamental, nsim=1000, output="network")

## 3.2.1. within-layer reciprocity

### 3.2.1.1. reciprocity in financial support layer

##### salaried
(emp.cash.sal.recip <-summary(cash_network~mutual:nodefactor("salaried"), interact.dependent="silent"))
sim.cash.sal.recip <-summary(sim.model~L(~mutual:nodefactor("salaried"),~cash), interact.dependent="silent")

##### status
(emp.cash.stat.recip <-summary(cash_network~mutual:nodefactor("status"), interact.dependent="silent"))
sim.cash.stat.recip <-summary(sim.model~L(~mutual:nodefactor("status"),~cash), interact.dependent="silent")

##### land
(emp.cash.land.recip <-summary(cash_network~mutual:nodefactor("land"), interact.dependent="silent"))
sim.cash.land.recip <-summary(sim.model~L(~mutual:nodefactor("land"),~cash), interact.dependent="silent")

##### calculating MC p value
pval1 <- apply(sim.cash.sal.recip <= emp.cash.sal.recip,2, mean)
pval1.top <- apply(sim.cash.sal.recip >= emp.cash.sal.recip,2, mean)
pval1 <- cbind(emp.cash.sal.recip,apply(sim.cash.sal.recip, 2,min), apply(sim.cash.sal.recip, 2,mean),
               apply(sim.cash.sal.recip, 2,max), pmin(1,2*pmin(pval1,pval1.top)))

pval2 <- apply(sim.cash.stat.recip <= emp.cash.stat.recip,2, mean)
pval2.top <- apply(sim.cash.stat.recip >= emp.cash.stat.recip,2, mean)
pval2 <- cbind(emp.cash.stat.recip,apply(sim.cash.stat.recip, 2,min), apply(sim.cash.stat.recip, 2,mean),
               apply(sim.cash.stat.recip, 2,max), pmin(1,2*pmin(pval2,pval2.top)))

pval3 <- apply(sim.cash.land.recip <= emp.cash.land.recip,2, mean)
pval3.top <- apply(sim.cash.land.recip >= emp.cash.land.recip,2, mean)
pval3 <- cbind(emp.cash.land.recip,apply(sim.cash.land.recip, 2,min), apply(sim.cash.land.recip, 2,mean),
               apply(sim.cash.land.recip, 2,max), pmin(1,2*pmin(pval3,pval3.top)))

pval <- rbind(pval3, pval1, pval2)

dimnames(pval)[[1]] <- c("land", "salaried", "status")
dimnames(pval)[[2]] <- c("obs","min","mean","max","MC p-value")
pval

##### plotting simulation results (nsim=1000)

cols <- c("grey","skyblue")

par(mfrow = c(3, 1), mar=c(5,3,3,3))

range <- min(emp.cash.land.recip,min(sim.cash.land.recip)):max(emp.cash.land.recip,max(sim.cash.land.recip))
freq <- tabulate(sim.cash.land.recip)
freq <- freq[range]
ci <- quantile(sim.cash.land.recip, probs = c(0.025, 0.975))
pos <- range >= min(ci) & range <= max(ci)

x<-barplot(freq, col = cols[pos + 1], border = cols[pos + 1], ylim=range(pretty(freq)),
           cex.axis = 2)
abline(v = x[emp.cash.land.recip-min(range)+1], col = "navy", lwd = 4, lty = 1)

range <- min(emp.cash.sal.recip,min(sim.cash.sal.recip)):max(emp.cash.sal.recip,max(sim.cash.sal.recip))
freq <- tabulate(sim.cash.sal.recip)
freq <- freq[range]
ci <- quantile(sim.cash.sal.recip, probs = c(0.025, 0.975))
pos <- range >= min(ci) & range <= max(ci)

x<-barplot(freq, col = cols[pos + 1], border = cols[pos + 1], ylim=range(pretty(freq)),
           cex.axis = 2)
abline(v = x[emp.cash.sal.recip-min(range)+1], col = "navy", lwd = 4, lty = 1)

range <- min(emp.cash.stat.recip,min(sim.cash.stat.recip)):max(emp.cash.stat.recip,max(sim.cash.stat.recip))
freq <- tabulate(sim.cash.stat.recip)
freq <- freq[range]
ci <- quantile(sim.cash.stat.recip, probs = c(0.025, 0.975))
pos <- range >= min(ci) & range <= max(ci)

x<-barplot(freq, col = cols[pos + 1], border = cols[pos + 1], ylim=range(pretty(freq)),
           cex.axis = 2)
abline(v = x[emp.cash.stat.recip-min(range)+1], col = "navy", lwd = 4, lty = 1)

### 3.2.1.2. reciprocity in material support layer

##### salaried
(emp.material.sal.recip <-summary(material_network~mutual:nodefactor("salaried"), interact.dependent="silent"))
sim.material.sal.recip <-summary(sim.model~L(~mutual:nodefactor("salaried"),~material), interact.dependent="silent")

##### status
(emp.material.stat.recip <-summary(material_network~mutual:nodefactor("status"), interact.dependent="silent"))
sim.material.stat.recip <-summary(sim.model~L(~mutual:nodefactor("status"),~material), interact.dependent="silent")

##### land
(emp.material.land.recip <-summary(material_network~mutual:nodefactor("land"), interact.dependent="silent"))
sim.material.land.recip <-summary(sim.model~L(~mutual:nodefactor("land"),~material), interact.dependent="silent")

##### calculating MC p value
pval1 <- apply(sim.material.sal.recip <= emp.material.sal.recip,2, mean)
pval1.top <- apply(sim.material.sal.recip >= emp.material.sal.recip,2, mean)
pval1 <- cbind(emp.material.sal.recip,apply(sim.material.sal.recip, 2,min), apply(sim.material.sal.recip, 2,mean),
               apply(sim.material.sal.recip, 2,max), pmin(1,2*pmin(pval1,pval1.top)))

pval2 <- apply(sim.material.stat.recip <= emp.material.stat.recip,2, mean)
pval2.top <- apply(sim.material.stat.recip >= emp.material.stat.recip,2, mean)
pval2 <- cbind(emp.material.stat.recip,apply(sim.material.stat.recip, 2,min), apply(sim.material.stat.recip, 2,mean),
               apply(sim.material.stat.recip, 2,max), pmin(1,2*pmin(pval2,pval2.top)))

pval3 <- apply(sim.material.land.recip <= emp.material.land.recip,2, mean)
pval3.top <- apply(sim.material.land.recip >= emp.material.land.recip,2, mean)
pval3 <- cbind(emp.material.land.recip,apply(sim.material.land.recip, 2,min), apply(sim.material.land.recip, 2,mean),
               apply(sim.material.land.recip, 2,max), pmin(1,2*pmin(pval3,pval3.top)))

pval <- rbind(pval3, pval1, pval2)

dimnames(pval)[[1]] <- c("land", "salaried", "status")
dimnames(pval)[[2]] <- c("obs","min","mean","max","MC p-value")
pval

##### plotting simulation results (nsim=1000)

par(mfrow = c(3, 1), mar=c(5,1,1,1))

range <- min(emp.material.land.recip,min(sim.material.land.recip)):max(emp.material.land.recip,max(sim.material.land.recip))
freq <- tabulate(sim.material.land.recip)
freq <- freq[range]
ci <- quantile(sim.material.land.recip, probs = c(0.025, 0.975))
pos <- range >= min(ci) & range <= max(ci)

x<-barplot(freq, col = cols[pos + 1], border = cols[pos + 1], ylim=range(pretty(freq)),
           cex.axis = 2)
abline(v = x[emp.material.land.recip-min(range)+1], col = "navy", lwd = 4, lty = 1)

range <- min(emp.material.sal.recip,min(sim.material.sal.recip)):max(emp.material.sal.recip,max(sim.material.sal.recip))
freq <- tabulate(sim.material.sal.recip)
freq <- freq[range]
ci <- quantile(sim.material.sal.recip, probs = c(0.025, 0.975))
pos <- range >= min(ci) & range <= max(ci)

x<-barplot(freq, col = cols[pos + 1], border = cols[pos + 1], ylim=range(pretty(freq)),
           cex.axis = 2)
abline(v = x[emp.material.sal.recip-min(range)+1], col = "navy", lwd = 4, lty = 1)

range <- min(emp.material.stat.recip,min(sim.material.stat.recip)):max(emp.material.stat.recip,max(sim.material.stat.recip))
freq <- tabulate(sim.material.stat.recip)
freq <- freq[range]
ci <- quantile(sim.material.stat.recip, probs = c(0.025, 0.975))
pos <- range >= min(ci) & range <= max(ci)

x<-barplot(freq, col = cols[pos + 1], border = cols[pos + 1], ylim=range(pretty(freq)),
           cex.axis = 2)
abline(v = x[emp.material.stat.recip-min(range)+1], col = "navy", lwd = 4, lty = 1)

### 3.2.1.3. reciprocity in labor support layer

##### salaried
(emp.help.sal.recip <-summary(help_network~mutual:nodefactor("new_salaried"), interact.dependent="silent"))
sim.help.sal.recip <-summary(sim.model~L(~mutual:nodefactor("new_salaried"),~help), interact.dependent="silent")

##### status
(emp.help.stat.recip <-summary(help_network~mutual:nodefactor("new_status"), interact.dependent="silent"))
sim.help.stat.recip <-summary(sim.model~L(~mutual:nodefactor("new_status"),~help), interact.dependent="silent")

##### land
(emp.help.land.recip <-summary(help_network~mutual:nodefactor("new_land2"), interact.dependent="silent"))
sim.help.land.recip <-summary(sim.model~L(~mutual:nodefactor("new_land2"),~help), interact.dependent="silent")


##### calculating MC p value
pval1 <- apply(sim.help.sal.recip <= emp.help.sal.recip,2, mean)
pval1.top <- apply(sim.help.sal.recip >= emp.help.sal.recip,2, mean)
pval1 <- cbind(emp.help.sal.recip,apply(sim.help.sal.recip, 2,min), apply(sim.help.sal.recip, 2,mean),
               apply(sim.help.sal.recip, 2,max), pmin(1,2*pmin(pval1,pval1.top)))

pval2 <- apply(sim.help.stat.recip <= emp.help.stat.recip,2, mean)
pval2.top <- apply(sim.help.stat.recip >= emp.help.stat.recip,2, mean)
pval2 <- cbind(emp.help.stat.recip,apply(sim.help.stat.recip, 2,min), apply(sim.help.stat.recip, 2,mean),
               apply(sim.help.stat.recip, 2,max), pmin(1,2*pmin(pval2,pval2.top)))

pval3 <- apply(sim.help.land.recip <= emp.help.land.recip,2, mean)
pval3.top <- apply(sim.help.land.recip >= emp.help.land.recip,2, mean)
pval3 <- cbind(emp.help.land.recip,apply(sim.help.land.recip, 2,min), apply(sim.help.land.recip, 2,mean),
               apply(sim.help.land.recip, 2,max), pmin(1,2*pmin(pval3,pval3.top)))

pval <- rbind(pval3, pval1, pval2)

dimnames(pval)[[1]] <- c("land", "salaried", "status")
dimnames(pval)[[2]] <- c("obs","min","mean","max","MC p-value")
pval

##### plotting simulation results (nsim=1000)

par(mfrow = c(3, 1), mar=c(5,1,1,1))

range <- min(emp.help.land.recip,min(sim.help.land.recip)):max(emp.help.land.recip,max(sim.help.land.recip))
freq <- tabulate(sim.help.land.recip+1)
ci <- quantile(sim.help.land.recip, probs = c(0.025, 0.975))
pos <- range >= min(ci) & range <= max(ci)

x<-barplot(freq, col = cols[pos + 1], border = cols[pos + 1], ylim=range(pretty(freq)),
           cex.axis = 2)
abline(v = x[emp.help.land.recip-min(range)+1], col = "navy", lwd = 4, lty = 1)

range <- min(emp.help.sal.recip,min(sim.help.sal.recip)):max(emp.help.sal.recip,max(sim.help.sal.recip))
freq <- tabulate(sim.help.sal.recip+1)
ci <- quantile(sim.help.sal.recip, probs = c(0.025, 0.975))
pos <- range >= min(ci) & range <= max(ci)

x<-barplot(freq, col = cols[pos + 1], border = cols[pos + 1], ylim=range(pretty(freq)),
           cex.axis = 2)
abline(v = x[emp.help.sal.recip-min(range)+1], col = "navy", lwd = 4, lty = 1)

range <- min(emp.help.stat.recip,min(sim.help.stat.recip)):max(emp.help.stat.recip,max(sim.help.stat.recip))
freq <- tabulate(sim.help.stat.recip+1)
ci <- quantile(sim.help.stat.recip, probs = c(0.025, 0.975))
pos <- range >= min(ci) & range <= max(ci)

x<-barplot(freq, col = cols[pos + 1], border = cols[pos + 1], ylim=range(pretty(freq)),
           cex.axis = 2)
abline(v = x[emp.help.stat.recip-min(range)+1], col = "navy", lwd = 4, lty = 1)

## 3.2.2. cross-layer reciprocity

### 3.2.2.1. financial-material support cross-layer reciprocity

##### salaried
(emp.cash.material.sal <- sum(summary(multi_network ~ L(~nodemix("salaried", levels=TRUE, levels2=c(2,4)), ~ cash & t(material)))))
sim.cash.material.sal <- rowSums(summary(sim.model ~ L(~nodemix("salaried", levels=TRUE, levels2=c(2,4)), ~ cash & t(material))))

##### status
(emp.cash.material.stat <- sum(summary(multi_network ~ L(~nodemix("status", levels=TRUE, levels2=c(2,4)), ~ cash & t(material)))))
sim.cash.material.stat <- rowSums(summary(sim.model ~ L(~nodemix("status", levels=TRUE, levels2=c(2,4)), ~ cash & t(material))))

##### land
(emp.cash.material.land <- sum(summary(multi_network ~ L(~nodemix("land", levels=TRUE, levels2=c(2,4)), ~ cash & t(material)))))
sim.cash.material.land <- rowSums(summary(sim.model ~ L(~nodemix("land", levels=TRUE, levels2=c(2,4)), ~ cash & t(material))))

##### calculating MC p value
sim.cash.material.sal <- as.matrix(sim.cash.material.sal)
sim.cash.material.stat <- as.matrix(sim.cash.material.stat)
sim.cash.material.land <- as.matrix(sim.cash.material.land)

pval1 <- apply(sim.cash.material.sal < emp.cash.material.sal,2, mean)
pval1.top <- apply(sim.cash.material.sal > emp.cash.material.sal,2, mean)
pval1 <- cbind(emp.cash.material.sal,apply(sim.cash.material.sal, 2,min), apply(sim.cash.material.sal, 2,mean),
               apply(sim.cash.material.sal, 2,max), pmin(1,2*pmin(pval1,pval1.top)))

pval2 <- apply(sim.cash.material.stat < emp.cash.material.stat,2, mean)
pval2.top <- apply(sim.cash.material.stat > emp.cash.material.stat,2, mean)
pval2 <- cbind(emp.cash.material.stat,apply(sim.cash.material.stat, 2,min), apply(sim.cash.material.stat, 2,mean),
               apply(sim.cash.material.stat, 2,max), pmin(1,2*pmin(pval2,pval2.top)))

pval3 <- apply(sim.cash.material.land < emp.cash.material.land,2, mean)
pval3.top <- apply(sim.cash.material.land > emp.cash.material.land,2, mean)
pval3 <- cbind(emp.cash.material.land,apply(sim.cash.material.land, 2,min), apply(sim.cash.material.land, 2,mean),
               apply(sim.cash.material.land, 2,max), pmin(1,2*pmin(pval3,pval3.top)))

pval <- rbind(pval3, pval1, pval2)

dimnames(pval)[[1]] <- c("land", "salaried", "status")
dimnames(pval)[[2]] <- c("obs","min","mean","max","MC p-value")
pval

##### plotting simulation results (nsim=1000)

par(mfrow = c(3, 1), mar=c(5,1,1,1))

range <- min(emp.cash.material.land, min(sim.cash.material.land)):max(emp.cash.material.land,max(sim.cash.material.land))
freq <- tabulate(sim.cash.material.land)
freq <- freq[range]
ci <- quantile(sim.cash.material.land, probs = c(0.025, 0.975))
pos <- range > min(ci) & range < max(ci)

x<-barplot(freq, col = cols[pos + 1], border = cols[pos + 1], ylim=range(pretty(freq)),
           cex.axis = 2)
abline(v = x[emp.cash.material.land-min(range)+1], col = "navy", lwd = 4, lty = 1)

range <- min(emp.cash.material.sal, min(sim.cash.material.sal)):max(emp.cash.material.sal,max(sim.cash.material.sal))
freq <- tabulate(sim.cash.material.sal)
freq <- freq[range]
ci <- quantile(sim.cash.material.sal, probs = c(0.025, 0.975))
pos <- range > min(ci) & range < max(ci)

x<-barplot(freq, col = cols[pos + 1], border = cols[pos + 1], ylim=range(pretty(freq)),
           cex.axis = 2)
abline(v = x[emp.cash.material.sal-min(range)+1], col = "navy", lwd = 4, lty = 1)

range <- min(emp.cash.material.stat, min(sim.cash.material.stat)):max(emp.cash.material.stat,max(sim.cash.material.stat))
freq <- tabulate(sim.cash.material.stat)
freq <- freq[range]
ci <- quantile(sim.cash.material.stat, probs = c(0.025, 0.975))
pos <- range > min(ci) & range < max(ci)

x<-barplot(freq, col = cols[pos + 1], border = cols[pos + 1], ylim=range(pretty(freq)),
           cex.axis = 2)
abline(v = x[emp.cash.material.stat-min(range)+1], col = "navy", lwd = 4, lty = 1)

### 3.2.2.2. financial-labor support cross-layer reciprocity 

##### salaried
(emp.cash.help.sal <- sum(summary(multi_network ~ L(~nodemix("salaried", levels=TRUE, levels2=c(2,4)), ~ cash & t(help)))))
sim.cash.help.sal <- rowSums(summary(sim.model ~ L(~nodemix("salaried", levels=TRUE, levels2=c(2,4)), ~ cash & t(help))))

##### status
(emp.cash.help.stat <- sum(summary(multi_network ~ L(~nodemix("status", levels=TRUE, levels2=c(2,4)), ~ cash & t(help)))))
sim.cash.help.stat <- rowSums(summary(sim.model ~ L(~nodemix("status", levels=TRUE, levels2=c(2,4)), ~ cash & t(help))))

##### land
(emp.cash.help.land <- sum(summary(multi_network ~ L(~nodemix("land", levels=TRUE, levels2=c(2,4)), ~ cash & t(help)))))
sim.cash.help.land <- rowSums(summary(sim.model ~ L(~nodemix("land", levels=TRUE, levels2=c(2,4)), ~ cash & t(help))))

##### calculating MC p value
sim.cash.help.sal <- as.matrix(sim.cash.help.sal)
sim.cash.help.stat <- as.matrix(sim.cash.help.stat)
sim.cash.help.land <- as.matrix(sim.cash.help.land)

pval1 <- apply(sim.cash.help.sal < emp.cash.help.sal,2, mean)
pval1.top <- apply(sim.cash.help.sal > emp.cash.help.sal,2, mean)
pval1 <- cbind(emp.cash.help.sal,apply(sim.cash.help.sal, 2,min), apply(sim.cash.help.sal, 2,mean),
               apply(sim.cash.help.sal, 2,max), pmin(1,2*pmin(pval1,pval1.top)))

sum(sim.cash.help.sal > emp.cash.help.sal)

pval2 <- apply(sim.cash.help.stat < emp.cash.help.stat,2, mean)
pval2.top <- apply(sim.cash.help.stat > emp.cash.help.stat,2, mean)
pval2 <- cbind(emp.cash.help.stat,apply(sim.cash.help.stat, 2,min), apply(sim.cash.help.stat, 2,mean),
               apply(sim.cash.help.stat, 2,max), pmin(1,2*pmin(pval2,pval2.top)))

pval3 <- apply(sim.cash.help.land < emp.cash.help.land,2, mean)
pval3.top <- apply(sim.cash.help.land > emp.cash.help.land,2, mean)
pval3 <- cbind(emp.cash.help.land,apply(sim.cash.help.land, 2,min), apply(sim.cash.help.land, 2,mean),
               apply(sim.cash.help.land, 2,max), pmin(1,2*pmin(pval3,pval3.top)))

pval <- rbind(pval3, pval1, pval2)

dimnames(pval)[[1]] <- c("land", "salaried", "status")
dimnames(pval)[[2]] <- c("obs","min","mean","max","MC p-value")
pval

##### plotting simulation results (nsim=1000)

par(mfrow = c(3, 1), mar=c(5,3,1,1))

range <- min(emp.cash.help.land, min(sim.cash.help.land)):max(emp.cash.help.land,max(sim.cash.help.land))
freq <- tabulate(sim.cash.help.land)
freq <- freq[range]
ci <- quantile(sim.cash.help.land, probs = c(0.025, 0.975))
pos <- range > min(ci) & range < max(ci)

x<-barplot(freq, col = cols[pos + 1], border = cols[pos + 1], ylim=range(pretty(freq)),
           cex.axis = 2)
abline(v = x[emp.cash.help.land-min(range)+1], col = "navy", lwd = 4, lty = 1)

range <- min(emp.cash.help.sal, min(sim.cash.help.sal)):max(emp.cash.help.sal,max(sim.cash.help.sal))
freq <- tabulate(sim.cash.help.sal)
freq <- freq[range]
ci <- quantile(sim.cash.help.sal, probs = c(0.025, 0.975))
pos <- range > min(ci) & range < max(ci)

x<-barplot(freq, col = cols[pos + 1], border = cols[pos + 1], ylim=range(pretty(freq)),
           cex.axis = 2)
abline(v = x[emp.cash.help.sal-min(range)+1], col = "navy", lwd = 4, lty = 1)

range <- min(emp.cash.help.stat, min(sim.cash.help.stat)):max(emp.cash.help.stat,max(sim.cash.help.stat))
freq <- tabulate(sim.cash.help.stat)
freq <- freq[range]
ci <- quantile(sim.cash.help.stat, probs = c(0.025, 0.975))
pos <- range > min(ci) & range < max(ci)

x<-barplot(freq, col = cols[pos + 1], border = cols[pos + 1], ylim=range(pretty(freq)),
           cex.axis = 2)
abline(v = x[emp.cash.help.stat-min(range)+1], col = "navy", lwd = 4, lty = 1)

### 3.2.2.3. material-labor support cross-layer reciprocity

##### salaried
(emp.material.help.sal <- sum(summary(multi_network ~ L(~nodemix("salaried", levels=TRUE, levels2=c(2,4)), ~ material & t(help)))))
sim.material.help.sal <- rowSums(summary(sim.model ~ L(~nodemix("salaried", levels=TRUE, levels2=c(2,4)), ~ material & t(help))))

##### status
(emp.material.help.stat <- sum(summary(multi_network ~ L(~nodemix("status", levels=TRUE, levels2=c(2,4)), ~ material & t(help)))))
sim.material.help.stat <- rowSums(summary(sim.model ~ L(~nodemix("status", levels=TRUE, levels2=c(2,4)), ~ material & t(help))))

##### land
(emp.material.help.land <- sum(summary(multi_network ~ L(~nodemix("land", levels=TRUE, levels2=c(2,4)), ~ material & t(help)))))
sim.material.help.land <- rowSums(summary(sim.model ~ L(~nodemix("land", levels=TRUE, levels2=c(2,4)), ~ material & t(help))))

##### calculating MC p value
sim.material.help.sal <- as.matrix(sim.material.help.sal)
sim.material.help.stat <- as.matrix(sim.material.help.stat)
sim.material.help.land <- as.matrix(sim.material.help.land)

pval1 <- apply(sim.material.help.sal < emp.material.help.sal,2, mean)
pval1.top <- apply(sim.material.help.sal > emp.material.help.sal,2, mean)
pval1 <- cbind(emp.material.help.sal,apply(sim.material.help.sal, 2,min), apply(sim.material.help.sal, 2,mean),
               apply(sim.material.help.sal, 2,max), pmin(1,2*pmin(pval1,pval1.top)))

pval2 <- apply(sim.material.help.stat < emp.material.help.stat,2, mean)
pval2.top <- apply(sim.material.help.stat > emp.material.help.stat,2, mean)
pval2 <- cbind(emp.material.help.stat,apply(sim.material.help.stat, 2,min), apply(sim.material.help.stat, 2,mean),
               apply(sim.material.help.stat, 2,max), pmin(1,2*pmin(pval2,pval2.top)))

pval3 <- apply(sim.material.help.land < emp.material.help.land,2, mean)
pval3.top <- apply(sim.material.help.land > emp.material.help.land,2, mean)
pval3 <- cbind(emp.material.help.land,apply(sim.material.help.land, 2,min), apply(sim.material.help.land, 2,mean),
               apply(sim.material.help.land, 2,max), pmin(1,2*pmin(pval3,pval3.top)))

pval <- rbind(pval3, pval1, pval2)

dimnames(pval)[[1]] <- c("land", "salaried", "status")
dimnames(pval)[[2]] <- c("obs","min","mean","max","MC p-value")
pval

##### plotting simulation results (nsim=1000)

par(mfrow = c(3, 1), mar=c(5,3,1,1))

range <- min(emp.material.help.land, min(sim.material.help.land)):max(emp.material.help.land,max(sim.material.help.land))
freq <- tabulate(sim.material.help.land)
freq <- freq[range]
ci <- quantile(sim.material.help.land, probs = c(0.025, 0.975))
pos <- range > min(ci) & range < max(ci)

x<-barplot(freq, col = cols[pos + 1], border = cols[pos + 1], ylim=range(pretty(freq)),
           cex.axis = 2)
abline(v = x[emp.material.help.land-min(range)+1], col = "navy", lwd = 4, lty = 1)

range <- min(emp.material.help.sal, min(sim.material.help.sal)):max(emp.material.help.sal,max(sim.material.help.sal))
freq <- tabulate(sim.material.help.sal)
freq <- freq[range]
ci <- quantile(sim.material.help.sal, probs = c(0.025, 0.975))
pos <- range > min(ci) & range < max(ci)

x<-barplot(freq, col = cols[pos + 1], border = cols[pos + 1], ylim=range(pretty(freq)),
           cex.axis = 2)
abline(v = x[emp.material.help.sal-min(range)+1], col = "navy", lwd = 4, lty = 1)

range <- min(emp.material.help.stat, min(sim.material.help.stat)):max(emp.material.help.stat,max(sim.material.help.stat))
freq <- tabulate(sim.material.help.stat)
freq <- freq[range]
ci <- quantile(sim.material.help.stat, probs = c(0.025, 0.975))
pos <- range > min(ci) & range < max(ci)

x<-barplot(freq, col = cols[pos + 1], border = cols[pos + 1], ylim=range(pretty(freq)),
           cex.axis = 2)
abline(v = x[emp.material.help.stat-min(range)+1], col = "navy", lwd = 4, lty = 1)

### 3.2.2.4. material-financial support cross-layer reciprocity (in reverse to 3.2.2.1)

##### salaried
(emp.material.cash.sal <- sum(summary(multi_network ~ L(~nodemix("salaried", levels=TRUE, levels2=c(3,4)), ~ cash & t(material)))))
sim.material.cash.sal <- rowSums(summary(sim.model ~ L(~nodemix("salaried", levels=TRUE, levels2=c(3,4)), ~ cash & t(material))))

##### status
(emp.material.cash.stat <- sum(summary(multi_network ~ L(~nodemix("status", levels=TRUE, levels2=c(3,4)), ~ cash & t(material)))))
sim.material.cash.stat <- rowSums(summary(sim.model ~ L(~nodemix("status", levels=TRUE, levels2=c(3,4)), ~ cash & t(material))))

##### land
(emp.material.cash.land <- sum(summary(multi_network ~ L(~nodemix("land", levels=TRUE, levels2=c(3,4)), ~ cash & t(material)))))
sim.material.cash.land <- rowSums(summary(sim.model ~ L(~nodemix("land", levels=TRUE, levels2=c(3,4)), ~ cash & t(material))))

##### calculating MC p value
sim.material.cash.sal <- as.matrix(sim.material.cash.sal)
sim.material.cash.stat <- as.matrix(sim.material.cash.stat)
sim.material.cash.land <- as.matrix(sim.material.cash.land)

pval1 <- apply(sim.material.cash.sal < emp.material.cash.sal,2, mean)
pval1.top <- apply(sim.material.cash.sal > emp.material.cash.sal,2, mean)
pval1 <- cbind(emp.material.cash.sal,apply(sim.material.cash.sal, 2,min), apply(sim.material.cash.sal, 2,mean),
               apply(sim.material.cash.sal, 2,max), pmin(1,2*pmin(pval1,pval1.top)))

pval2 <- apply(sim.material.cash.stat < emp.material.cash.stat,2, mean)
pval2.top <- apply(sim.material.cash.stat > emp.material.cash.stat,2, mean)
pval2 <- cbind(emp.material.cash.stat,apply(sim.material.cash.stat, 2,min), apply(sim.material.cash.stat, 2,mean),
               apply(sim.material.cash.stat, 2,max), pmin(1,2*pmin(pval2,pval2.top)))

pval3 <- apply(sim.material.cash.land < emp.material.cash.land,2, mean)
pval3.top <- apply(sim.material.cash.land > emp.material.cash.land,2, mean)
pval3 <- cbind(emp.material.cash.land,apply(sim.material.cash.land, 2,min), apply(sim.material.cash.land, 2,mean),
               apply(sim.material.cash.land, 2,max), pmin(1,2*pmin(pval3,pval3.top)))

pval <- rbind(pval3, pval1, pval2)

dimnames(pval)[[1]] <- c("land", "salaried", "status")
dimnames(pval)[[2]] <- c("obs","min","mean","max","MC p-value")
pval

##### plotting simulation results (nsim=1000)

par(mfrow = c(3, 1), mar=c(5,3,1,1))

range <- min(emp.material.cash.land, min(sim.material.cash.land)):max(emp.material.cash.land,max(sim.material.cash.land))
freq <- tabulate(sim.material.cash.land)
freq <- freq[range]
ci <- quantile(sim.material.cash.land, probs = c(0.025, 0.975))
pos <- range > min(ci) & range < max(ci)

x<-barplot(freq, col = cols[pos + 1], border = cols[pos + 1], ylim=range(pretty(freq)),
           cex.axis = 2)
abline(v = x[emp.material.cash.land-min(range)+1], col = "navy", lwd = 4, lty = 1)

range <- min(emp.material.cash.sal, min(sim.material.cash.sal)):max(emp.material.cash.sal,max(sim.material.cash.sal))
freq <- tabulate(sim.material.cash.sal)
freq <- freq[range]
ci <- quantile(sim.material.cash.sal, probs = c(0.025, 0.975))
pos <- range > min(ci) & range < max(ci)

x<-barplot(freq, col = cols[pos + 1], border = cols[pos + 1], ylim=range(pretty(freq)),
           cex.axis = 2)
abline(v = x[emp.material.cash.sal-min(range)+1], col = "navy", lwd = 4, lty = 1)

range <- min(emp.material.cash.stat, min(sim.material.cash.stat)):max(emp.material.cash.stat,max(sim.material.cash.stat))
freq <- tabulate(sim.material.cash.stat)
freq <- freq[range]
ci <- quantile(sim.material.cash.stat, probs = c(0.025, 0.975))
pos <- range > min(ci) & range < max(ci)

x<-barplot(freq, col = cols[pos + 1], border = cols[pos + 1], ylim=range(pretty(freq)),
           cex.axis = 2)
abline(v = x[emp.material.cash.stat-min(range)+1], col = "navy", lwd = 4, lty = 1)


### 3.2.2.5. labor-financial support cross-layer reciprocity (in reverse to 3.2.2.2)

##### salaried
(emp.help.cash.sal <- sum(summary(multi_network ~ L(~nodemix("salaried", levels=TRUE, levels2=c(3,4)), ~ cash & t(help)))))
sim.help.cash.sal <- rowSums(summary(sim.model ~ L(~nodemix("salaried", levels=TRUE, levels2=c(3,4)), ~ cash & t(help))))

##### status
(emp.help.cash.stat <- sum(summary(multi_network ~ L(~nodemix("status", levels=TRUE, levels2=c(3,4)), ~ cash & t(help)))))
sim.help.cash.stat <- rowSums(summary(sim.model ~ L(~nodemix("status", levels=TRUE, levels2=c(3,4)), ~ cash & t(help))))

##### land
(emp.help.cash.land <- sum(summary(multi_network ~ L(~nodemix("land", levels=TRUE, levels2=c(3,4)), ~ cash & t(help)))))
sim.help.cash.land <- rowSums(summary(sim.model ~ L(~nodemix("land", levels=TRUE, levels2=c(3,4)), ~ cash & t(help))))

##### calculating MC p value
sim.help.cash.sal <- as.matrix(sim.help.cash.sal)
sim.help.cash.stat <- as.matrix(sim.help.cash.stat)
sim.help.cash.land <- as.matrix(sim.help.cash.land)

pval1 <- apply(sim.help.cash.sal < emp.help.cash.sal,2, mean)
pval1.top <- apply(sim.help.cash.sal > emp.help.cash.sal,2, mean)
pval1 <- cbind(emp.help.cash.sal,apply(sim.help.cash.sal, 2,min), apply(sim.help.cash.sal, 2,mean),
               apply(sim.help.cash.sal, 2,max), pmin(1,2*pmin(pval1,pval1.top)))

pval2 <- apply(sim.help.cash.stat < emp.help.cash.stat,2, mean)
pval2.top <- apply(sim.help.cash.stat > emp.help.cash.stat,2, mean)
pval2 <- cbind(emp.help.cash.stat,apply(sim.help.cash.stat, 2,min), apply(sim.help.cash.stat, 2,mean),
               apply(sim.help.cash.stat, 2,max), pmin(1,2*pmin(pval2,pval2.top)))

pval3 <- apply(sim.help.cash.land < emp.help.cash.land,2, mean)
pval3.top <- apply(sim.help.cash.land > emp.help.cash.land,2, mean)
pval3 <- cbind(emp.help.cash.land,apply(sim.help.cash.land, 2,min), apply(sim.help.cash.land, 2,mean),
               apply(sim.help.cash.land, 2,max), pmin(1,2*pmin(pval3,pval3.top)))

pval <- rbind(pval3, pval1, pval2)

dimnames(pval)[[1]] <- c("land", "salaried", "status")
dimnames(pval)[[2]] <- c("obs","min","mean","max","MC p-value")
pval

##### plotting simulation results (nsim=1000)

par(mfrow = c(3, 1), mar=c(5,3,1,1))

range <- min(emp.help.cash.land, min(sim.help.cash.land)):max(emp.help.cash.land,max(sim.help.cash.land))
freq <- tabulate(sim.help.cash.land)
freq <- freq[range]
ci <- quantile(sim.help.cash.land, probs = c(0.025, 0.975))
pos <- range > min(ci) & range < max(ci)

x<-barplot(freq, col = cols[pos + 1], border = cols[pos + 1], ylim=range(pretty(freq)),
           cex.axis = 2)
abline(v = x[emp.help.cash.land-min(range)+1], col = "navy", lwd = 4, lty = 1)

range <- min(emp.help.cash.sal, min(sim.help.cash.sal)):max(emp.help.cash.sal,max(sim.help.cash.sal))
freq <- tabulate(sim.help.cash.sal)
freq <- freq[range]
ci <- quantile(sim.help.cash.sal, probs = c(0.025, 0.975))
pos <- range > min(ci) & range < max(ci)

x<-barplot(freq, col = cols[pos + 1], border = cols[pos + 1], ylim=range(pretty(freq)),
           cex.axis = 2)
abline(v = x[emp.help.cash.sal-min(range)+1], col = "navy", lwd = 4, lty = 1)

range <- min(emp.help.cash.stat, min(sim.help.cash.stat)):max(emp.help.cash.stat,max(sim.help.cash.stat))
freq <- tabulate(sim.help.cash.stat)
freq <- freq[range]
ci <- quantile(sim.help.cash.stat, probs = c(0.025, 0.975))
pos <- range > min(ci) & range < max(ci)

x<-barplot(freq, col = cols[pos + 1], border = cols[pos + 1], ylim=range(pretty(freq)),
           cex.axis = 2)
abline(v = x[emp.help.cash.stat-min(range)+1], col = "navy", lwd = 4, lty = 1)

### 3.2.2.6. labor-material support cross-layer reciprocity (in reverse to 3.2.2.3)

##### salaried
(emp.help.material.sal <- sum(summary(multi_network ~ L(~nodemix("salaried", levels=TRUE, levels2=c(3,4)), ~ material & t(help)))))
sim.help.material.sal <- rowSums(summary(sim.model ~ L(~nodemix("salaried", levels=TRUE, levels2=c(3,4)), ~ material & t(help))))

##### status
(emp.help.material.stat <- sum(summary(multi_network ~ L(~nodemix("status", levels=TRUE, levels2=c(3,4)), ~ material & t(help)))))
sim.help.material.stat <- rowSums(summary(sim.model ~ L(~nodemix("status", levels=TRUE, levels2=c(3,4)), ~ material & t(help))))

##### land
(emp.help.material.land <- sum(summary(multi_network ~ L(~nodemix("land", levels=TRUE, levels2=c(3,4)), ~ material & t(help)))))
sim.help.material.land <- rowSums(summary(sim.model ~ L(~nodemix("land", levels=TRUE, levels2=c(3,4)), ~ material & t(help))))

##### calculating MC p value
sim.help.material.sal <- as.matrix(sim.help.material.sal)
sim.help.material.stat <- as.matrix(sim.help.material.stat)
sim.help.material.land <- as.matrix(sim.help.material.land)

pval1 <- apply(sim.help.material.sal < emp.help.material.sal,2, mean)
pval1.top <- apply(sim.help.material.sal > emp.help.material.sal,2, mean)
pval1 <- cbind(emp.help.material.sal,apply(sim.help.material.sal, 2,min), apply(sim.help.material.sal, 2,mean),
               apply(sim.help.material.sal, 2,max), pmin(1,2*pmin(pval1,pval1.top)))

pval2 <- apply(sim.help.material.stat < emp.help.material.stat,2, mean)
pval2.top <- apply(sim.help.material.stat > emp.help.material.stat,2, mean)
pval2 <- cbind(emp.help.material.stat,apply(sim.help.material.stat, 2,min), apply(sim.help.material.stat, 2,mean),
               apply(sim.help.material.stat, 2,max), pmin(1,2*pmin(pval2,pval2.top)))

pval3 <- apply(sim.help.material.land < emp.help.material.land,2, mean)
pval3.top <- apply(sim.help.material.land > emp.help.material.land,2, mean)
pval3 <- cbind(emp.help.material.land,apply(sim.help.material.land, 2,min), apply(sim.help.material.land, 2,mean),
               apply(sim.help.material.land, 2,max), pmin(1,2*pmin(pval3,pval3.top)))

pval <- rbind(pval3, pval1, pval2)

dimnames(pval)[[1]] <- c("land", "salaried", "status")
dimnames(pval)[[2]] <- c("obs","min","mean","max","MC p-value")
pval

##### plotting simulation results (nsim=1000)

par(mfrow = c(3, 1), mar=c(5,3,1,1))

range <- min(emp.help.material.land, min(sim.help.material.land)):max(emp.help.material.land,max(sim.help.material.land))
freq <- tabulate(sim.help.material.land)
freq <- freq[range]
ci <- quantile(sim.help.material.land, probs = c(0.025, 0.975))
pos <- range > min(ci) & range < max(ci)

x<-barplot(freq, col = cols[pos + 1], border = cols[pos + 1], ylim=range(pretty(freq)),
           cex.axis = 2)
abline(v = x[emp.help.material.land-min(range)+1], col = "navy", lwd = 4, lty = 1)

range <- min(emp.help.material.sal, min(sim.help.material.sal)):max(emp.help.material.sal,max(sim.help.material.sal))
freq <- tabulate(sim.help.material.sal)
freq <- freq[range]
ci <- quantile(sim.help.material.sal, probs = c(0.025, 0.975))
pos <- range > min(ci) & range < max(ci)

x<-barplot(freq, col = cols[pos + 1], border = cols[pos + 1], ylim=range(pretty(freq)),
           cex.axis = 2)
abline(v = x[emp.help.material.sal-min(range)+1], col = "navy", lwd = 4, lty = 1)

range <- min(emp.help.material.stat, min(sim.help.material.stat)):max(emp.help.material.stat,max(sim.help.material.stat))
freq <- tabulate(sim.help.material.stat)
freq <- freq[range]
ci <- quantile(sim.help.material.stat, probs = c(0.025, 0.975))
pos <- range > min(ci) & range < max(ci)

x<-barplot(freq, col = cols[pos + 1], border = cols[pos + 1], ylim=range(pretty(freq)),
           cex.axis = 2)
abline(v = x[emp.help.material.stat-min(range)+1], col = "navy", lwd = 4, lty = 1)

