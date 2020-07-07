library(bayou)

tree <- read.tree("tree.tree")
mydata3 <- read.csv("mydata3.csv")
mydata3[,1]  
setdiff(tree$tip.label, mydata3[,1])
setdiff(mydata3[,1],tree$tip.label)

#pca1
dat.pca1 <- mydata3[,2]
names(dat.pca1) <- as.character(mydata3[,1])

priorOU.pca1 <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dk="cdpois", dtheta="dnorm"), 
                            param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                       dk=list(lambda=4, kmax=8), dsb=list(bmax=1, prob=1), 
                                       dtheta=list(mean=mean(dat.pca1), sd=0.5)),plot.prior = FALSE)





startpars.pca1 <- priorSim(priorOU.pca1, tree, plot=FALSE, cex=0.5)$pars[[1]]
priorOU.pca1(startpars.pca1)




set.seed(1)
#run1
mcmcOU.pca1_1 <- bayou.makeMCMC(tree, dat.pca1, SE=0.5, prior=priorOU.pca1, 
                               new.dir=TRUE, outname="./modelOU_rpca1_1", plot.freq=NULL) # Set up the MCMC
mcmcOU.pca1_1$run(2500000) # Run the MCMC
chainOU.pca1_1 <- mcmcOU.pca1_1$load()
summary.pca1_1 <- summary(chainOU.pca1_1)
save(chainOU.pca1_1, file = "mcmcOU.pca1_1.Rdata")

#run2
mcmcOU.pca1_2 <- bayou.makeMCMC(tree, dat.pca1, SE=0.5, prior=priorOU.pca1, 
                               new.dir=TRUE, outname="./modelOU_rpca1_2", plot.freq=NULL) # Set up the MCMC
mcmcOU.pca1_2$run(2500000) # Run the MCMC
chainOU.pca1_2 <- mcmcOU.pca1_2$load()
summary.pca1_2 <- summary(chainOU.pca1_2)
save(chainOU.pca1_2, file = "mcmcOU.pca1_2.Rdata")

#run3
mcmcOU.pca1_3 <- bayou.makeMCMC(tree, dat.pca1, SE=0.5, prior=priorOU.pca1, 
                               new.dir=TRUE, outname="./modelOU_rpca1_3", plot.freq=NULL) # Set up the MCMC
mcmcOU.pca1_3$run(2500000) # Run the MCMC
chainOU.pca1_3 <- mcmcOU.pca1_3$load()
summary.pca1_3 <- summary(chainOU.pca1_3)
save(chainOU.pca1_3, file = "mcmcOU.pca1_3.Rdata")


list_chains.pca1 <- list(chainOU.pca1_1,chainOU.pca1_2,chainOU.pca1_3)

combine.pca1 <- combine.chains(list_chains.pca1, burnin.prop = 0.3 )

summary.pca1 <- summary(combine.pca1)

plotSimmap.mcmc(combine.pca1, burnin = 0.3, pp.cutoff = 0.3, cex = 0.5)
plotBranchHeatMap(tree, combine.pca1, "theta", burnin = 0.3, pal = rainbow, cex = 0.5)
phenogram.density(tree, dat.pca1, burnin = 0.3, combine.pca1, pp.cutoff = 0.3)
par(mfrow=c(3, 5) )
plot(combine.pca1, auto.layout=FALSE)


