library(bayou)

tree <- read.tree("tree.tree")
mydata3 <- read.csv("mydata3.csv")
mydata3[,1]  
setdiff(tree$tip.label, mydata3[,1])
setdiff(mydata3[,1],tree$tip.label)

#pca3
dat.pca3 <- mydata3[,4]
names(dat.pca3) <-as.character(mydata3[,1])

priorOU.pca3 <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dk="cdpois", dtheta="dnorm"), 
                            param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                       dk=list(lambda=4, kmax=8), dsb=list(bmax=1, prob=1), 
                                       dtheta=list(mean=mean(dat.pca3), sd=0.5)),plot.prior = FALSE)





startpars.pca3 <- priorSim(priorOU.pca3, tree, plot=FALSE, cex=0.5)$pars[[1]]
priorOU.pca3(startpars.pca3)




set.seed(1)
#run1
mcmcOU.pca3_1 <- bayou.makeMCMC(tree, dat.pca3, SE=0.5, prior=priorOU.pca3, 
                               new.dir=TRUE, outname="./modelOU_rpca3_1", plot.freq=NULL) # Set up the MCMC
mcmcOU.pca3_1$run(2500000) # Run the MCMC
chainOU.pca3_1 <- mcmcOU.pca3_1$load()
summary.pca3_1 <- summary(chainOU.pca3_1)
save(chainOU.pca3_1, file = "mcmcOU.pca3_1.Rdata")

#run2
mcmcOU.pca3_2 <- bayou.makeMCMC(tree, dat.pca3, SE=0.5, prior=priorOU.pca3, 
                               new.dir=TRUE, outname="./modelOU_rpca3_2", plot.freq=NULL) # Set up the MCMC
mcmcOU.pca3_2$run(2500000) # Run the MCMC
chainOU.pca3_2 <- mcmcOU.pca3_2$load()
summary.pca3_2 <- summary(chainOU.pca3_2)
save(chainOU.pca3_2, file = "mcmcOU.pca3_2.Rdata")

#run3
mcmcOU.pca3_3 <- bayou.makeMCMC(tree, dat.pca3, SE=0.5, prior=priorOU.pca3, 
                               new.dir=TRUE, outname="./modelOU_rpca3_3", plot.freq=NULL) # Set up the MCMC
mcmcOU.pca3_3$run(2500000) # Run the MCMC
chainOU.pca3_3 <- mcmcOU.pca3_3$load()
summary.pca3_3 <- summary(chainOU.pca3_3)
save(chainOU.pca3_3, file = "mcmcOU.pca3_3.Rdata")


list_chains.pca3 <- list(chainOU.pca3_1,chainOU.pca3_2,chainOU.pca3_3)

combine.pca3 <- combine.chains(list_chains.pca3, burnin.prop = 0.3 )

summary.pca3 <- summary(combine.pca3)

plotSimmap.mcmc(combine.pca3, burnin = 0.3, pp.cutoff = 0.3, cex = 0.5)
plotBranchHeatMap(tree, combine.pca3, "theta", burnin = 0.3, pal = rainbow, cex = 0.5)
phenogram.density(tree, dat.pca3, burnin = 0.3, combine.pca3, pp.cutoff = 0.3)
par(mfrow=c(3, 5) )
plot(combine.pca3, auto.layout=FALSE)


