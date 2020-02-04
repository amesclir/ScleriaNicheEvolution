library(bayou)

tree <- read.tree("tree.tree")
mydata <- read.csv("mydata.csv")
mydata2 <- read.csv("mydata2.csv")
mydata[,1]  
setdiff(tree$tip.label, mydata[,1])
setdiff(mydata[,1],tree$tip.label)

#bio15
dat.bio15 <- mydata[,16]
names(dat.bio15) <- mydata[,1]
se.bio15 <- sqrt(mydata2[,16])
names(se.bio15) <- mydata2[,1]

priorOU.bio15 <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dk="cdpois", dtheta="dnorm"), 
                            param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                       dk=list(lambda=4, kmax=8), dsb=list(bmax=1, prob=1), 
                                       dtheta=list(mean=mean(dat.bio15), sd=1.5*sd(dat.bio15))),plot.prior = FALSE)





startpars.bio15 <- priorSim(priorOU.bio15, tree, plot=FALSE, cex=0.5)$pars[[1]]
priorOU.bio15(startpars.bio15)




set.seed(1)
#run1
mcmcOU.bio15_1 <- bayou.makeMCMC(tree, dat.bio15, SE=se.bio15, prior=priorOU.bio15, 
                               new.dir=TRUE, outname="./modelOU_rbio15_1", plot.freq=NULL) # Set up the MCMC
mcmcOU.bio15_1$run(2500000) # Run the MCMC
chainOU.bio15_1 <- mcmcOU.bio15_1$load()
summary.bio15_1 <- summary(chainOU.bio15_1)
save(chainOU.bio15_1, file = "mcmcOU.bio15_1.Rdata")

#run2
mcmcOU.bio15_2 <- bayou.makeMCMC(tree, dat.bio15, SE=se.bio15, prior=priorOU.bio15, 
                               new.dir=TRUE, outname="./modelOU_rbio15_2", plot.freq=NULL) # Set up the MCMC
mcmcOU.bio15_2$run(2500000) # Run the MCMC
chainOU.bio15_2 <- mcmcOU.bio15_2$load()
summary.bio15_2 <- summary(chainOU.bio15_2)
save(chainOU.bio15_2, file = "mcmcOU.bio15_2.Rdata")

#run3
mcmcOU.bio15_3 <- bayou.makeMCMC(tree, dat.bio15, SE=se.bio15, prior=priorOU.bio15, 
                               new.dir=TRUE, outname="./modelOU_rbio15_3", plot.freq=NULL) # Set up the MCMC
mcmcOU.bio15_3$run(2500000) # Run the MCMC
chainOU.bio15_3 <- mcmcOU.bio15_3$load()
summary.bio15_3 <- summary(chainOU.bio15_3)
save(chainOU.bio15_3, file = "mcmcOU.bio15_3.Rdata")


list_chains.bio15 <- list(chainOU.bio15_1,chainOU.bio15_2,chainOU.bio15_3)

combine.bio15 <- combine.chains(list_chains.bio15, burnin.prop = 0.3 )

summary.bio15 <- summary(combine.bio15)

plotSimmap.mcmc(combine.bio15, burnin = 0.3, pp.cutoff = 0.3, cex = 0.5)
plotBranchHeatMap(tree, combine.bio15, "theta", burnin = 0.3, pal = rainbow, cex = 0.5)
phenogram.density(tree, dat.bio15, burnin = 0.3, combine.bio15, pp.cutoff = 0.3)
par(mfrow=c(3, 5) )
plot(combine.bio15, auto.layout=FALSE)


