####figure6 network and properteis
library(igraph)
library(dplyr)
library(Hmisc)
library(qgraph)

col_g <- "#C1C1C1"
cols <- c("#DEB99B" ,"#5ECC6D", "#5DAFD9", "#7ED1E4", "#EA9527", "#F16E1D" ,"#6E4821", "#A4B423",
          "#C094DF" ,"#DC95D8" ,"#326530", "#50C0C9", "#67C021" ,"#DC69AF", "#8C384F", "#30455C", "#F96C72","#5ED2BF")

cormatrix.input  <- matrix(0, 364, 364)
##download from the MENA
cormatrix.input[row(cormatrix.input) >= col(cormatrix.input)] <- scan("soil_Correlation.txt")

cormatrix <- t(cormatrix.input)

for (i in 1:nrow(cormatrix)){
  for (j in 1:ncol(cormatrix)){
    if (i>j){cormatrix[i,j]<-cormatrix[j,i]}
  }
}
###the data file
name <- read.table("soil_species.txt", sep='\t',h=F,row.names=1,check=F,comment='')

name1 <- row.names(name)

rownames(cormatrix) <- name1

colnames(cormatrix) <- name1

###obtian from the MENA
cormatrix[abs(cormatrix)< 0.93] <- 0

g <- graph.adjacency(cormatrix, weighted = TRUE, mode = 'undirected')

g <- simplify(g)

g <- delete.vertices(g, which(degree(g)==0) )


E(g)$correlation <- E(g)$weight

E(g)$weight <- abs(E(g)$weight)

set.seed(007)

V(g)$modularity <- membership(cluster_fast_greedy(g))

V(g)$label <- V(g)$name

V(g)$label <- NA

modu_sort <- V(g)$modularity %>% table() %>% sort(decreasing = T)

top_num <- 10

modu_name <- names(modu_sort[1:10])

modu_cols <- cols[1:length(modu_name)]

names(modu_cols) <- modu_name

V(g)$color <- V(g)$modularity

V(g)$color[!(V(g)$color %in% modu_name)] <- col_g

V(g)$color[(V(g)$color %in% modu_name)] <- modu_cols[match(V(g)$color[(V(g)$color %in% modu_name)],modu_name)]

V(g)$frame.color <- V(g)$color


E(g)$color <- col_g

for ( i in modu_name){
  col_edge <- cols[which(modu_name==i)]
  otu_same_modu <-V(g)$name[which(V(g)$modularity==i)]
  E(g)$color[(data.frame(as_edgelist(g))$X1 %in% otu_same_modu)&(data.frame(as_edgelist(g))$X2 %in% otu_same_modu)] <- col_edge
}

e <- get.edgelist(g,names=FALSE)

r <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(g),
                                       area=8*(vcount(g)^2),repulse.rad=(vcount(g)^3.1))

par(font.main=4)

plot(g,layout=r, edge.color = E(g)$color,vertex.size=2)

title(main = paste0('Nodes=',length(V(g)$name),', ','Edges=',nrow(data.frame(as_edgelist(g)))))

####robustness
data <- read.csv("soil.csv", header = T, row.names = 1)
counts<-rowSums(data>0)

otutab<-data[counts>=33,]

comm<-t(otutab)
###calculate the relative abundance
sp.ra<-colMeans(comm)/15010

cormatrix=matrix(0,ncol(comm),ncol(comm))

for (i in 1:ncol(comm)){
  for (j in i:ncol(comm)){
    speciesi<-sapply(1:nrow(comm),function(k){
      ifelse(comm[k,i]>0,comm[k,i],ifelse(comm[k,j]>0,0.01,NA))
    })
    speciesj<-sapply(1:nrow(comm),function(k){
      ifelse(comm[k,j]>0,comm[k,j],ifelse(comm[k,i]>0,0.01,NA))
    })
    corij<-cor(log(speciesi)[!is.na(speciesi)],log(speciesj)[!is.na(speciesj)])
    
    cormatrix[i,j]<-cormatrix[j,i]<-corij
    
  }}

row.names(cormatrix)<-colnames(cormatrix)<-colnames(comm) # if processed using MENAP, OTU order should match in the original OTU table and the correlation matrix downloaded from MENAP.

cormatrix2<-cormatrix*(abs(cormatrix)>=0.91)  #only keep links above the cutoff point
cormatrix2[is.na(cormatrix2)]<-0
diag(cormatrix2)<-0    #no links for self-self    
sum(abs(cormatrix2)>0)/2  #this should be the number of links. 
sum(colSums(abs(cormatrix2))>0)  # node number: number of species with at least one linkage with others.

network.raw<-cormatrix2[colSums(abs(cormatrix2))>0,colSums(abs(cormatrix2))>0]
sp.ra2<-sp.ra[colSums(abs(cormatrix2))>0]
sum(row.names(network.raw)==names(sp.ra2))  #check if matched



rand.remov.once<-function(netRaw, rm.percent, sp.ra, abundance.weighted=T){
  id.rm<-sample(1:nrow(netRaw), round(nrow(netRaw)*rm.percent))
  net.Raw=netRaw #don't want change netRaw
  net.Raw[id.rm,]=0;  net.Raw[,id.rm]=0;   ##remove all the links to these species
  if (abundance.weighted){
    net.stength= net.Raw*sp.ra
  } else {
    net.stength= net.Raw
  }
  
  sp.meanInteration<-colMeans(net.stength)
  
  id.rm2<- which(sp.meanInteration<=0)  ##remove species have negative interaction or no interaction with others
  remain.percent<-(nrow(netRaw)-length(id.rm2))/nrow(netRaw)
  #for simplicity, I only consider the immediate effects of removing the
  #'id.rm' species; not consider the sequential effects of extinction of
  # the 'id.rm2' species.
  
  #you can write out the network pruned
  #  net.Raw[id.rm2,]=0;  net.Raw[,id.rm2]=0;
  # write.csv( net.Raw,"network pruned.csv")
  
  remain.percent
}

rm.p.list=seq(0.05,0.2,by=0.05)
rmsimu<-function(netRaw, rm.p.list, sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov.once(netRaw=netRaw, rm.percent=x, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
    })
    remain.mean=mean(remains)
    remain.sd=sd(remains)
    remain.se=sd(remains)/(nperm^0.5)
    result<-c(remain.mean,remain.sd,remain.se)
    names(result)<-c("remain.mean","remain.sd","remain.se")
    result
  }))
}

Weighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
Unweighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=F,nperm=100)

dat1<-data.frame(Proportion.removed=rep(seq(0.05,1,by=0.05),2),rbind(Weighted.simu,Unweighted.simu),
                 weighted=rep(c("weighted","unweighted"),each=20),
                 year=rep("ys",40),treat=rep("hpla",40))

currentdat<-dat1

write.csv(currentdat,"robustness.soil.csv")

######cohesion 
####################create necessary functions######################

#find the number of zeroes in a vector
zero <- function(vec){
  num.zero <- length(which(vec == 0))
  return(num.zero)
}

#create function that averages only negative values in a vector
neg.mean <- function(vector){
  neg.vals <- vector[which(vector < 0)]
  n.mean <- mean(neg.vals)
  if(length(neg.vals) == 0) n.mean <- 0
  return(n.mean)
}

#create function that averages only positive values in a vector
pos.mean <- function(vector){
  pos.vals <- vector[which(vector > 0)]
  p.mean <- mean(pos.vals)
  if(length(pos.vals) == 0) p.mean <- 0
  return(p.mean)
}

###################################################################
###################################################################
### Workflow options ####
###################################################################
###################################################################

## Choose a persistence cutoff (min. fraction of taxon presence) for retaining taxa in the analysis
pers.cutoff <- 0.09
## Decide the number of iterations to run for each taxon. (>= 200 is recommended)
# Larger values of iter mean the script takes longer to run
iter <- 600
## Decide whether to use taxon/column shuffle (tax.shuffle = T) or row shuffle algorithm (tax.shuffle = F)
tax.shuffle <- T
## Option to input your own correlation table
# Note that your correlation table MUST have the same number of taxa as the abundance table. There should be no empty (all zero) taxon vectors in the abundance table. 
# Even if you input your own correlation table, the persistence cutoff will be applied
use.custom.cors <- F

soil <- read.table("soil.txt", header = T,row.names=1,sep="\t")
b <- t(soil)
c <- as.matrix(b)
c <- c[rowSums(c) > 0, colSums(c) > 0]
rowsums.orig <- rowSums(c)
zero.cutoff <- ceiling(pers.cutoff * dim(c)[1])
d <- c[ , apply(c, 2, zero) < (dim(c)[1]-zero.cutoff) ]
d <- d[rowSums(d) > 0, ]
rel.d <- d / rowsums.orig
cor.mat.true <- cor(rel.d)

med.tax.cors <- vector()
if(use.custom.cors == F) {
  if(tax.shuffle) {
    for(which.taxon in 1:dim(rel.d)[2]){
      perm.cor.vec.mat <- vector()
      
      for(i in 1:iter){
        perm.rel.d <- matrix(numeric(0), dim(rel.d)[1], dim(rel.d)[2])
        rownames(perm.rel.d) <- rownames(rel.d)
        colnames(perm.rel.d) <- colnames(rel.d)
        
        #For each otu
        for(j in 1:dim(rel.d)[2]){ 
          # Replace the original taxon vector with a permuted taxon vector
          perm.rel.d[, j ] <- sample(rel.d[ ,j ]) 
        }
        
        # Do not randomize focal column 
        perm.rel.d[, which.taxon] <- rel.d[ , which.taxon]
        
        # Calculate correlation matrix of permuted matrix
        cor.mat.null <- cor(perm.rel.d)
        
        # For each iteration, save the vector of null matrix correlations between focal taxon and other taxa
        perm.cor.vec.mat <- cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])
        
      }
      # Save the median correlations between the focal taxon and all other taxa  
      med.tax.cors <- cbind(med.tax.cors, apply(perm.cor.vec.mat, 1, median))
      
      # For large datasets, this can be helpful to know how long this loop will run
      if(which.taxon %% 20 == 0){print(which.taxon)}
    }
  } else {
    for(which.taxon in 1:dim(rel.d)[2]){
      
      #create vector to hold correlations from every permutation for each single otu
      ## perm.cor.vec.mat stands for permuted correlations vector matrix
      perm.cor.vec.mat <- vector()
      
      for(i in 1:iter){
        #Create duplicate matrix to shuffle abundances
        perm.rel.d <- rel.d 
        
        #For each taxon
        for(j in 1:dim(rel.d)[1]){ 
          which.replace <- which(rel.d[j, ] > 0 ) 
          # if the focal taxon is greater than zero, take it out of the replacement vector, so the focal abundance stays the same
          which.replace.nonfocal <- which.replace[!(which.replace %in% which.taxon)]
          
          #Replace the original taxon vector with a vector where the values greater than 0 have been randomly permuted 
          perm.rel.d[j, which.replace.nonfocal] <- sample(rel.d[ j, which.replace.nonfocal]) 
        }
        
        # Calculate correlation matrix of permuted matrix
        cor.mat.null <- cor(perm.rel.d)
        
        # For each iteration, save the vector of null matrix correlations between focal taxon and other taxa
        perm.cor.vec.mat <- cbind(perm.cor.vec.mat, cor.mat.null[, which.taxon])
        
      }
      # Save the median correlations between the focal taxon and all other taxa  
      med.tax.cors <- cbind(med.tax.cors, apply(perm.cor.vec.mat, 1, median))
      
      # For large datasets, this can be helpful to know how long this loop will run
      if(which.taxon %% 20 == 0){print(which.taxon)}
    }
  }
}

ifelse(use.custom.cors == T, {
  obs.exp.cors.mat <- custom.cor.mat.sub}, {
    obs.exp.cors.mat <- cor.mat.true - med.tax.cors
  }
)

diag(obs.exp.cors.mat) <- 0

connectedness.pos <- apply(obs.exp.cors.mat, 2, pos.mean)
connectedness.neg <- apply(obs.exp.cors.mat, 2, neg.mean)

cohesion.pos <- rel.d %*% connectedness.pos
cohesion.neg <- rel.d %*% connectedness.neg

output <- list(connectedness.neg, connectedness.pos, cohesion.neg, cohesion.pos)
names(output) <- c("Negative Connectedness", "Positive Connectedness", "Negative Cohesion", "Positive Cohesion")

soil.co <- output
write.csv(cohesion.neg, file = "soil.co.neg.csv")
write.csv(cohesion.pos, file = "soil.co.pos.csv")
