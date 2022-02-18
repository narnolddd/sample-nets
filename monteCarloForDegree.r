#
# The original graph degree is what we want to approximate using Bayes'
#
# we observe sampled graph and use Monte Carlo with the assumpiton that we know the sample technique
# and the original number of nodes.
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

# Using a power law network

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
for(mno in 1:length(prob) ){
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

  # approximate the total number of links in g from gR, simple approximation 
  numLinksOri <- as.integer( sum(degR)/p )

  # approximate the degree sequenc of the original nodes from the degree of the observed graph
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

# ++++++++++++++++++++++++++++++++++++++++++++++++++++
# Monte Carlo: to approximate degG from degApp
# starting from gApp, remove links, evaluate the
# degree sequence and compare it with degR
# vectorizing the function
# two possible cases the average degree from the sample reconstructed network is equal to the degree of gR
# other case the degree of the sample reconstructed network is equal to the degree of gR
# the initial value resOld

# (code vectorised as runs a bit faster, instead of two loops)
  lenDegApp <- length(degApp)
  lstPred <- rep(0,lenDegApp) # this is the degree sequence when sampling gApp with probability p
  for(i in 1:lenDegApp){  
    seqJ <- seq(0,degApp[i])
    lstPred[i] <- sum(ourBinomial(degApp[i],seqJ,p)*seqJ)
  }
  
  # starting value
  resOld <- mean( (lstPred-degR)^2 )

#ts <- Sys.time() #timing the process
#lstConv <- NULL # too check how fast is the convergence

  # the randomisation and reconnections
  for(abc in 1:15000){
    posNodes <- sample(lenDegApp,2,replace=F) # takes two nodes at random
    
    #if a node 1 is of degree one do not use as we are going to remove a link from it and we don't want degree zero nodes
    while(degApp[posNodes[1]] <= 1){  
      posNodes <- sample(lenDegApp,2,replace=F)
    }
    degApp[posNodes[1]] <- degApp[posNodes[1]]-1 # decrease degree of node 1
    degApp[posNodes[2]] <- degApp[posNodes[2]]+1 # incerase degree of node 2

    lstPred <- rep(0,lenDegApp)
    for(i in 1:lenDegApp){  
      seqJ <- seq(0,degApp[i])
      lstPred[i] <- sum(ourBinomial(degApp[i],seqJ,p)*seqJ)
    }
    resNew <- mean( (lstPred-degR)^2 )
    if(resNew > resOld){ # go back to what we use to have reject step
      degApp[posNodes[1]] <- degApp[posNodes[1]]+1
      degApp[posNodes[2]] <- degApp[posNodes[2]]-1
    }else{ # accept step
      resOld <- resNew
#       lstConv <- cbind(lstConv,c(abc,resOld)) # for debugging
    }

  } # for abc
#  te <- Sys.time()
# +++++++++++++++++++ end Monte Carlo+++++++++++++++++++++++++++++++++

# mean square error of the apprximated degree with the original degree
  msqGGap <- mean( (degG - degApp)^2 )

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

  resExpe <- rbind(resExpe, c(p,msqGGap,msqGBayesAppro,msqGBayesExact) )
  resNumLinks <- rbind(resNumLinks, c(p, sum(degApp), sum(degBayesAppro), sum(degBayesExact)))
}

pdf("bayesDegApproERJan2022.pdf")
plot(resExpe[,1],log(resExpe[,4]),type = 'o',main="degAppr (blk exact) (blue app) (orange MC)") 
lines(resExpe[,1],log(resExpe[,2]),col='blue')
lines(resExpe[,1],log(resExpe[,3]),col='orange')
dev.off()

pdf("bayesLinkApproERJan2022.pdf")
plot(resNumLinks[,1],log(resNumLinks[,4]),type = 'o',main="linkAppr (blk exact) (blue app) (orange MC)") 
lines(resNumLinks[,1],log(resNumLinks[,2]),col='blue')
lines(resNumLinks[,1],log(resNumLinks[,3]),col='orange')
dev.off()

q("no")
