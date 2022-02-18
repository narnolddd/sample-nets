#
# Approximate the prior by
# 1 ranking the nodes in decreasing degree
# 2 use the basic approximation k_i ~ k'_i/p
# 3 redistribute the degree between the nodes such that there are k'_i=0
# Conserve the total number of degrees
#

library('igraph')

# ++++++++++++++++++++++++++++++++++++++++++++++++++++
# Binomial approximation of sampling the nodes' degree
# ++++++++++++++++++++++++++++++++++++++++++++++++++++
ourBinomial <- function(ts,te,q)
{
  tmp1 <- choose(ts,te)*q^te*(1-q)^(ts-te)
  return(tmp1)
}
# ++++++++++++++++++++++++++++++++++++++++++++++++++++

# Using a power law graph
# g is the original graph (BA m=5) and degG its degree sequence, numLinks number of edges in original graph
#g <- simplify(sample_pa(1000,m=5,directed=F))
g <- g<-erdos.renyi.game(1000, 1/20,directed=F)

degG <- degree(g)
numLinksG <- ecount(g)

# probability that a link is observed is p,
# prob is a list with probabilities to try 
prob <- seq(1,20)/21.0
resExpe <- NULL
resNumLinks <- NULL
for(mno in 1:20){#length(prob) ){
  p <- prob[mno] #choose a probability
  cat("doing p= ", p ,"\n")
  pC <- 1-p                  # probability that link survives
  deltaLinks <- pC*ecount(g) # number of links that are observed from the original graph

  # select the links to be removed,
  linksToRemove <- sample(numLinksG,deltaLinks,replace=F);

  # new graph gR with the removed links, This is the graph that we observed, from which we are going to
  # try to reproduce the degree sequence of g from the degree sequence degR of the gR
  # as we are constructing gR from g we inherit the same number of nodes, even if they are isolated nodes (degree=0)
  gR <- delete_edges(g,linksToRemove)
  degR <- degree(gR)
  lstGR <- as_edgelist(gR)

  # approximate the total number of links in g from gR, simple approximation 
  numLinksOri <- as.integer( sum(degR)/p )

  # approximate the degree sequence of the original nodes from the degree of the observed graph
  degApp <- as.integer(degR/p)

  # as the sum of degApp is not equal to the predicted number of links we adjust degApp
  numMissLinks <- numLinksOri - sum(degApp)
  while(numMissLinks > 0){ # we need to add links until there are no missing links
    ranNode <- sample(length(degApp),numMissLinks,replace=F) # select numMissLinks random nodes
    degApp[ranNode] <- degApp[ranNode] + 1 # increase their degree
    numMissLinks <- numLinksOri - sum(degApp)
  }
  numMissLinks <- numLinksOri - sum(degApp)
  while(numMissLinks < 0){ # we need to remove links until there are no missing links
    ranNode <- sample(length(degApp),numMissLinks,replace=F) # select numMissLinks random nodes
    ind <- which(degApp[ranNode] != 0) # only select nodes that do not have zero degree as we are going to remove links
    degApp[ranNode[ind]] <- degApp[ranNode] - 1 # decrease their degree
    numMissLinks <- numLinksOri - sum(degApp)
  }

  # +++++++++++++++++++++++++ part to rank the nodes and create the new degree sequence  +++++++
  # sort the degrees and create an index vector reference to degApp
  resSort <- sort.int(degApp,decreasing=TRUE,index.return=TRUE)
  # the invIndx is the index to return to the original labelling of the nodes
  invIndx <- sort.int(resSort$ix,decreasing=FALSE, index.return=TRUE)

  # make the adjustment using a copy of degApp
  # notice that the ranking is not necessary but I like it (probably for later on for rich-club)
  degTmp <- degApp[resSort$ix] 
  ind <- min(which(degTmp == 0))
  while(length(ind) != 0){
    degTmp[ind-1] <- degTmp[ind-1]-1
    degTmp[ind] <- degTmp[ind]+1
    minNod <- which(degTmp <= 0)
    if(length(minNod) != 0){
      ind <- min(minNod)
    }else{
      ind <- minNod
    }
  }
  # back to the original ranking
  degApp <- degTmp[invIndx$ix]
  

# as we have change degApp, evaluate degApp-degR to evaluate how many links are missing per node
  degExtra <- degApp - degR; # if all is correct degExtra >=0 but check just in case
  if(sum(degExtra)%%2 == 1){ #if this is odd then we cannot create a graph change it to even number
    degExtra[3] <- degExtra[3]+1;
  }
# construct a new graph with the extra links, we simplify to remove loops and multilinks and get edgelist
  gOld <- simplify(sample_degseq(degExtra,method=c("simple")))
  lstGN <- as_edgelist(gOld)
  
# the new graph is formed by the edge list of gR and the edge list of gN
  lstNew <- rbind(lstGR,lstGN)
  gApp <- simplify(graph_from_edgelist(lstNew,directed=F))
  degApp <- degree(gApp) # the new degree sequence

# ++++++++++++++++++ Bayes ++++++++++++++++++++++++++++++++++

# approximating the degree sequence using degApp as the prior
  priorDeg <- degApp/sum(degApp); # prior for degree sequence
  probDeno1 <- NULL
  degBayesAppro <- rep(0,lenDegApp)
  for(i in 1:lenDegApp){
    deno <- sum(ourBinomial(degApp[i:lenDegApp],degR[i],p)*priorDeg[i:lenDegApp])
    nume <- sum(ourBinomial(degApp[i:lenDegApp],degR[i],p)*priorDeg[i:lenDegApp]*degApp[i:lenDegApp])
    degBayesAppro[i] <- nume/deno
  }
  msqGBayesAppro <- mean( (degG-degBayesAppro)^2 )

# approximating the degree sequence using degG as the prior
  priorDeg <- degG/sum(degG); # prior for degree sequence
  probDeno1 <- NULL
  degBayesExact <- rep(0,lenDegApp)
  for(i in 1:lenDegApp){
    deno <- sum(ourBinomial(degApp[i:lenDegApp],degR[i],p)*priorDeg[i:lenDegApp])
    nume <- sum(ourBinomial(degApp[i:lenDegApp],degR[i],p)*priorDeg[i:lenDegApp]*degApp[i:lenDegApp])
    degBayesExact[i] <- nume/deno
  }
  msqGBayesExact <- mean( (degG-degBayesExact)^2 )

  msqAppro <- mean( (degG-degApp)^2 )
  resExpe <- rbind(resExpe, c(p, msqAppro, msqGBayesAppro, msqGBayesExact) )
  
} # for mno

plot(resExpe[,1],log(resExpe[,4]),type = 'o',main="degAppr (blk exact) (blue app) (orange MC)",ylim=c(-1,6)) 
lines(resExpe[,1],log(resExpe[,2]),col='blue')
lines(resExpe[,1],log(resExpe[,3]),col='orange')
#dev.off()

#pdf("reuseBayesLinkApproBAJan2022.pdf")
#plot(resNumLinks[,1],(resNumLinks[,2]/sum(degG)),type = 'o',main="linkAppr(blck app)(blue app)(orange MC)",ylim=c(0.5,1.5)) 
#lines(resNumLinks[,1],(resNumLinks[,3]/sum(degG)),col='blue')
#lines(resNumLinks[,1],(resNumLinks[,4]/sum(degG)),col='orange')
#dev.off()

#q("no")
