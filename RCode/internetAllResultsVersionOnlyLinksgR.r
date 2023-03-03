# from the networks generate by Naomi
# evaluate the approximation to the number of triangles 

library('igraph')
source("NullModelRichClub.r")
source("kPlusEvaluation.r")
source("degAndKPlusNewLabel.r")
source("kPlusForAppro.r")
source("probEnsemble.r")

ourBinomial <- function(ts,te,q)
{
  tmp1 <- choose(ts,te)*q^te*(1-q)^(ts-te)
  return(tmp1)
}

#
# properties of the real network
#
dat <- read.table("rutasInter.dat")
res <- cbind(dat$V1,dat$V2)
g <- graph_from_edgelist(res,directed = F)
numNodesOriginal <- vcount(g)
numNodes <- numNodesOriginal
numLinksOriginal <- ecount(g)
degG <- degree(g)
# triangles per node
numTriangles<-count_triangles(g)
# triangles per link
matA <- get.adjacency(g)
linkA <- as_edgelist(g)
triangle_matrixA <- as.matrix(Matrix::tcrossprod(matA))
#the links for graph g
triangleLink <- triangle_matrixA[linkA]
linkTria <- cbind(linkA,triangleLink)
totTria <- sum(triangleLink);

# check that 2T_i = \sum _j TiJ
checkNumTria <- rep(0,numNodes)
for(i in 1:(length(linkTria[,3]))){
  checkNumTria[linkTria[i,1]] <- checkNumTria[linkTria[i,1]]+linkTria[i,3]
  checkNumTria[linkTria[i,2]] <- checkNumTria[linkTria[i,2]]+linkTria[i,3]  
}
#the check is that checkNumTria/2 = numTriangles
write.table(as_edgelist(g),"InterNaomi/REAL",col.names=F,row.names=F);
# +++ the sampled graphs
lstResults <- NULL
probList <- rev( c(0.1,0.2,0.3,0.4,0.45,0.5,0.55,0.6,0.7,0.8,0.9) );

#for(ip in 1:(length(probList))){
for(ip in 6:(length(probList))){
prob <- probList[ip];
lstResAll <- NULL
lstResLinks <- NULL
lstTotTria <- NULL
fileLabMax <- 30

for(fileLab in 1:fileLabMax){
 cat("doing ", fileLab, "\n")
  p <- prob;
  pC <- 1-p
  deltaLinks <- pC*ecount(g)
  linksToRemove <- sample(ecount(g),deltaLinks,replace=F);
  gR <- delete_edges(g,linksToRemove)
  degR <- degree(gR)
  nameFile <- sprintf("InterNaomi/INTER%2.1f-%d",prob,fileLab);
  write.table(as_edgelist(gR),nameFile,col.names=F,row.names=F);

# check if g and gR have the same number of nodes, if not
  # to evaluate only the links in the gR graph
  matB <- get.adjacency(gR)
  linkB <- as_edgelist(gR)
  
# add some to gR so we can do the montecarlo
  diffNodes <- as.integer(vcount(g)-vcount(gR))
  if(diffNodes>0){
    gR <- add_vertices(gR,diffNodes)
    degR <- degree(gR)
  }

# +++++++++++++++++++++++++++++++++++++++++++++
# evaluate the degree distribution using Bayes
# +++++++++++++++++++++++++++++++++++++++++++++
  priorDeg <- degG/sum(degG); # prior for degree sequence
  probDeno1 <- NULL
  degBayesAppro <- rep(0,vcount(g))
  for(i in 1:20)  degBayesAppro[i] <- degR[i]/p
  for(i in 21:vcount(g)){
    indSel <- which(degG >= degR[i]) # new part
    deno <- sum(ourBinomial(degG[indSel],degR[i],p)*priorDeg[indSel])
    nume <- sum(ourBinomial(degG[indSel],degR[i],p)*priorDeg[indSel]*degG[indSel])
    degBayesAppro[i] <- nume/deno
#    if(nume == Inf | deno == Inf) degBayesAppro[i] = degR[i]/p;#as the binomial explodes
#    if(is.nan(nume) | is.nan(deno)) degBayesAppro[i] = degR[i]/p;#as the binomial explodes
  }
  degBayesAppro[1]  <- degR[1]/p
  
# convert it to integers
  degBayesInt <- as.integer((degBayesAppro/sum(degBayesAppro))*sum(degG));
  numLinksBayes <- sum(degBayesInt);

# new version May 18th
# missing number of links
  degSeq <- degBayesInt

  missingLinks = as.integer((sum(degR)/p-sum(degSeq)))
  cat("missingLinks ", missingLinks, "\n")
  while(missingLinks > 0){
#  if(missingLinks > 0){
    ind <- length(degSeq)
    if(missingLinks > ind) missingLinks <- ind
    newDeg <- sample(length(degSeq),missingLinks,replace=F)

# generate the ensemble

# add the missing links
    degSeq[newDeg] <- degSeq[newDeg]+1;
      missingLinks = as.integer((sum(degR)/p-sum(degSeq)))
      cat("missingLinksA ", missingLinks, "\n")	
}

#if(missingLinks < 0){
while(missingLinks <0){
 ind <- which(degSeq > 1)
  if(-missingLinks > length(ind)) missingLinks <- -length(ind)
 newDeg <- sample(length(ind),-missingLinks,replace=F)
 ind2 <- ind[newDeg]
 # generate the ensemble
#  degSeq <- degBayesInt
# remove the extra links
  degSeq[newDeg] <- degSeq[newDeg]-1;
        missingLinks = as.integer((sum(degR)/p-sum(degSeq)))
      cat("missingLinksB ", missingLinks, "\n")	
#stop()
}


if(sum(degSeq)%%2 == 1){ # we have an impossible to build graph, need to add a link
    degSeq[1] <- degSeq[1]+1
  }

# sort the degrees and keep track of the new labels
  degApproSorted <- sort.int(degSeq,decreasing=TRUE,index.return=TRUE);
  ind <- degApproSorted$ix
  degKplusGr <- degAndKPlusNewLabel(gR, degApproSorted)
  kPlusForApproSorted <- kPlusForAppro(degKplusGr,degApproSorted,numNodes);
  matProb <- probEnsemble(degApproSorted,kPlusForApproSorted)

## number of triangles per node
#  lstTrianEnsembleNew <- rep(0,numNodes)
#  for(i in 1:numNodes){
#    tmp <- 0
#    pij <- matProb[i,]
#    tmp1 <- 0
#    for(k in 1:numNodes){
#      pjk <- matProb[,k]
#      pki <- matProb[k,i]
#      tmp1 <- tmp1 + pjk*pki
#    }
#    tmp2 <- sum(tmp1*pij)/3
#    lstTrianEnsembleNew[i] <- tmp2
#  cat("i= ", i," tri= ",numTriangles[i]," appT= ", tmp2, "\n")
#  }
  
  #errorTriaNodes <- mean( (lstTrianEnsembleNew-numTriangles)^2)
  #sdTriaNodes <- sd( (lstTrianEnsembleNew-numTriangles)^2)
  # end count triangles per node
#  lstResAll <- cbind(lstResAll,errorTriaNodes)

# LINKS

  lstTriaLinks <- rep(0,numNodes*numNodes)
  lstLinksEntry <- NULL #rep(0,numNodes*numNodes)

#new method multiplying the matrix
  lstApproTria <- NULL
  for(a in 1:length(linkB[,1])){
    i <- linkB[a,1]
    j <- linkB[a,2]
    pij <- matProb[i,j]
    pik <- sum(matProb[j,]*matProb[,i])
    #old version links from the real network
    lstApproTria <- rbind(lstApproTria,c(triangle_matrixA[i,j],pij*pik))	
  }

  #old method multiplying the matrix
#  lstApproTria <- NULL
#  for(a in 1:length(linkTria[,1])){
#    i <- linkTria[a,1]
#    j <- linkTria[a,2]
#    pij <- matProb[i,j]
#    pik <- sum(matProb[j,]*matProb[,i])
#    lstApproTria <- rbind(lstApproTria,c(linkTria[a,3],pij*pik))
#  }

#   avTotTria <- mean( (2*sum(lstApproTria[,2])-totTria)^2)
#   lstTotTria <- cbind(lstTotTria,avTotTria)
#  lstResLinks <- cbind(lstResLinks, mean((lstApproTria[,1]-lstApproTria[,2])^2))
#cat("errorLinks ", mean((lstApproTria[,1]-lstApproTria[,2])^2),"avTotTria = ", avTotTria, "\n")
#  lstResLinks <- cbind(lstResLinks,linkApproTria[,2])
  # alternative method
#stop();

lstResLinks <- rbind(lstResLinks, c(mean((lstApproTria[,1]-lstApproTria[,2])^2), sd((lstApproTria[,1]-lstApproTria[,2])^2)))
cat("filelab ", fileLab, "numLinks ",sum(degSeq), " ori ", sum(degG),"\n");

} #for fileLaba

lstResults <- rbind(lstResults,c(p,mean(lstResLinks[,1]),mean(lstResLinks[,2])) )

#lstResults <- rbind(lstResults,c(p,mean(lstResLinks),sd(lstResLinks),sum(matProb)/2, numLinksOriginal,mean(lstTotTria),sd(lstTotTria)));
#cat("results ",p," ", mean(lstResLinks), " ", sd(lstResLinks), " ",sum(matProb)/2, " ",numLinksOriginal,"\n");
write.table(lstResults,"resultsINTEREnsembleTrianglesV4.dat",row.names=F,col.names=F)


}# for ip
write.table(lstResults,"resultsINTEREnsembleTrianglesV4.dat",row.names=F,col.names=F)

q("no")

stop(0)

#average by measurement that is av of fileLabMax then av of numnodes
avRes <- rep(0,numNodes)
sdRes <- rep(0,numNodes)
for(i in 1:numNodes){
  avRes[i] <- mean(lstResAll[i,])
  sdRes[i] <- sd(lstResAll[i,])
}
errorTriaNodesPer <- mean( (avRes-numTriangles)^2)
sdTriaNodesPer <- sd( (avRes-numTriangles)^2)

#cat("errorTriaPer ", errorTriaNodesPer, "\n");
#average by numNodes then by fileLabMax
avRes <- rep(0,fileLabMax)
sdRes <- rep(0,fileLabMax)
for(i in 1:fileLabMax){
  avRes[i] <- mean( (lstResAll[,i]-numTriangles)^2)
  sdRes[i] <- sd( (lstResAll[,i] -numTriangles)^2)
}

errorTriaNodes <- mean( avRes)
sdTriaNodes <- mean( sdRes) 

 avResLinks <- rep(0,fileLabMax)
 sdResLinks <- rep(0,fileLabMax)
 for(i in 1:fileLabMax){
 avResLinks[i] <- mean((lstResLinks[,i]-linkTria[,3])^2)
 sdResLinks[i] <- sd((lstResLinks[,i]-linkTria[,3])^2)
 }

#stop();
#plot(avRes,type='o')
#lines(sdRes)

#pdf("numTrianglesNode.pdf")
#plot(numTriangles[0:100],main=sprintf("p=%2.1f",prob))
#lines(avRes,col='red')
#lines(avRes+sdRes,col='blue')
#lines(avRes-sdRes,col='green')

# rescale the results
#r1 <- (avRes/sum(avRes))*sum(numTriangles)
#plot(numTriangles[0:100],main=sprintf("rescale p=%2.1f",prob))
#lines(r1,col='red');

# two ways to do the average ??
avResLinks <- rep(0,numLinksOriginal)
sdResLinks <- rep(0,numLinksOriginal)
for(i in 1:numLinksOriginal){
  avResLinks[i] <- mean(lstResLinks[i,])
  sdResLinks[i] <- sd(lstResLinks[i,])
}

errorTriaLinksPer <- mean((linkTria[,3]-avResLinks)^2)
sdTriaLinksPer <- sd((linkTria[,3]-avResLinks)^2)

 avResLinks <- rep(0,fileLabMax)
 sdResLinks <- rep(0,fileLabMax)
  totTria <- rep(0,49)
 sumTria <- sum(linkTria[,3])
 for(i in 1:fileLabMax){
 avResLinks[i] <- mean((lstResLinks[,i]-linkTria[,3])^2)
 sdResLinks[i] <- sd((lstResLinks[,i]-linkTria[,3])^2)
  totTria[i] <- (sum(lstResLinks[,i]) - sumTria)^2
 }
 errorTriaLinks <- mean(avResLinks)
 sdTriaLinks <- mean(sdResLinks)
  errorTotTria <- mean(totTria)
 sdTotTria <- sd(totTria)

cat("errorTriaNodes ", errorTriaNodesPer, " ",sdTriaNodesPer," ",errorTriaNodes," ",sdTriaNodes,"\n");
cat("errorTriaLinks ", errorTriaLinksPer," ", sdTriaLinksPer," ",errorTriaLinks," ",sdTriaLinks,"\n")
#stop();

#plot(linkTria[1:200,3],main=sprintf("T_ell p=%2.1f err=%4.1f",prob,errorTriaLinks))
#lines(avResLinks,col='red')
#lines(avResLinks+sdResLinks,col='blue')
#lines(avResLinks-sdResLinks,col='green')

# rescale the results
#r1 <- (avResLinks/sum(avResLinks))*sum(linkTria[,3])
#plot(linkTria[0:200,3],main=sprintf("rescale T_ell p=%2.1f",prob))
#lines(r1,col='red');

#dev.off()
lstResults <- rbind(lstResults, c(p,errorTriaNodesPer,sdTriaNodesPer,errorTriaNodes,sdTriaNodes,errorTriaLinksPer,sdTriaLinksPer,errorTriaLinks,sdTriaLinks,errorTotTria,sdTotTria));

#} #for ip

write.table(lstResults,"resultsINTEREnsembleTrianglesV2.dat",row.names=F,col.names=F)
q("no")




