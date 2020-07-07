library(bayou)

tree <- read.tree("tree.tree")
mydata3 <- read.csv("mydata3.csv")
mydata3[,1]  
setdiff(tree$tip.label, mydata3[,1])
setdiff(mydata3[,1],tree$tip.label)

#pca2
dat.pca2 <- mydata3[,3]
names(dat.pca2) <-as.character(mydata3[,1])

priorOU.pca2 <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dk="cdpois", dtheta="dnorm"), 
                            param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                       dk=list(lambda=4, kmax=8), dsb=list(bmax=1, prob=1), 
                                       dtheta=list(mean=mean(dat.pca2), sd=0.5)),plot.prior = FALSE)





startpars.pca2 <- priorSim(priorOU.pca2, tree, plot=FALSE, cex=0.5)$pars[[1]]
priorOU.pca2(startpars.pca2)




set.seed(1)
#run1
mcmcOU.pca2_1 <- bayou.makeMCMC(tree, dat.pca2, SE=0.5, prior=priorOU.pca2, 
                               new.dir=TRUE, outname="./modelOU_rpca2_1", plot.freq=NULL) # Set up the MCMC
mcmcOU.pca2_1$run(2500000) # Run the MCMC
chainOU.pca2_1 <- mcmcOU.pca2_1$load()
summary.pca2_1 <- summary(chainOU.pca2_1)
save(chainOU.pca2_1, file = "mcmcOU.pca2_1.Rdata")

#run2
mcmcOU.pca2_2 <- bayou.makeMCMC(tree, dat.pca2, SE=0.5, prior=priorOU.pca2, 
                               new.dir=TRUE, outname="./modelOU_rpca2_2", plot.freq=NULL) # Set up the MCMC
mcmcOU.pca2_2$run(2500000) # Run the MCMC
chainOU.pca2_2 <- mcmcOU.pca2_2$load()
summary.pca2_2 <- summary(chainOU.pca2_2)
save(chainOU.pca2_2, file = "mcmcOU.pca2_2.Rdata")

#run3
mcmcOU.pca2_3 <- bayou.makeMCMC(tree, dat.pca2, SE=0.5, prior=priorOU.pca2, 
                               new.dir=TRUE, outname="./modelOU_rpca2_3", plot.freq=NULL) # Set up the MCMC
mcmcOU.pca2_3$run(2500000) # Run the MCMC
chainOU.pca2_3 <- mcmcOU.pca2_3$load()
summary.pca2_3 <- summary(chainOU.pca2_3)
save(chainOU.pca2_3, file = "mcmcOU.pca2_3.Rdata")


list_chains.pca2 <- list(chainOU.pca2_1,chainOU.pca2_2,chainOU.pca2_3)

combine.pca2 <- combine.chains(list_chains.pca2, burnin.prop = 0.3 )

summary.pca2 <- summary(combine.pca2)

plotSimmap.mcmc(combine.pca2, burnin = 0.3, pp.cutoff = 0.3, cex = 0.5)
plotBranchHeatMap(tree, combine.pca2, "theta", burnin = 0.3, pal = rainbow, cex = 0.5)
phenogram.density(tree, dat.pca2, burnin = 0.3, combine.pca2, pp.cutoff = 0.3)
par(mfrow=c(3, 5) )
plot(combine.pca2, auto.layout=FALSE)


