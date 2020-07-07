library(bayou)

tree <- read.tree("tree.tree")
mydata <- read.csv("mydata.csv")
mydata2 <- read.csv("mydata2.csv")
mydata[,1]  
setdiff(tree$tip.label, mydata[,1])
setdiff(mydata[,1],tree$tip.label)

#bio12
dat.bio12 <- mydata[,13]
names(dat.bio12) <- mydata[,1]
se.bio12 <- sqrt(mydata2[,13])
names(se.bio12) <- mydata2[,1]

priorOU.bio12 <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dk="cdpois", dtheta="dnorm"), 
                            param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                       dk=list(lambda=4, kmax=8), dsb=list(bmax=1, prob=1), 
                                       dtheta=list(mean=mean(dat.bio12), sd=1.5*sd(dat.bio12))),plot.prior = FALSE)





startpars.bio12 <- priorSim(priorOU.bio12, tree, plot=FALSE, cex=0.5)$pars[[1]]
priorOU.bio12(startpars.bio12)




set.seed(1)
#run1
mcmcOU.bio12_1 <- bayou.makeMCMC(tree, dat.bio12, SE=se.bio12, prior=priorOU.bio12, 
                               new.dir=TRUE, outname="./modelOU_rbio12_1", plot.freq=NULL) # Set up the MCMC
mcmcOU.bio12_1$run(2500000) # Run the MCMC
chainOU.bio12_1 <- mcmcOU.bio12_1$load()
summary.bio12_1 <- summary(chainOU.bio12_1)
save(chainOU.bio12_1, file = "mcmcOU.bio12_1.Rdata")

#run2
mcmcOU.bio12_2 <- bayou.makeMCMC(tree, dat.bio12, SE=se.bio12, prior=priorOU.bio12, 
                               new.dir=TRUE, outname="./modelOU_rbio12_2", plot.freq=NULL) # Set up the MCMC
mcmcOU.bio12_2$run(2500000) # Run the MCMC
chainOU.bio12_2 <- mcmcOU.bio12_2$load()
summary.bio12_2 <- summary(chainOU.bio12_2)
save(chainOU.bio12_2, file = "mcmcOU.bio12_2.Rdata")

#run3
mcmcOU.bio12_3 <- bayou.makeMCMC(tree, dat.bio12, SE=se.bio12, prior=priorOU.bio12, 
                               new.dir=TRUE, outname="./modelOU_rbio12_3", plot.freq=NULL) # Set up the MCMC
mcmcOU.bio12_3$run(2500000) # Run the MCMC
chainOU.bio12_3 <- mcmcOU.bio12_3$load()
summary.bio12_3 <- summary(chainOU.bio12_3)
save(chainOU.bio12_3, file = "mcmcOU.bio12_3.Rdata")


list_chains.bio12 <- list(chainOU.bio12_1,chainOU.bio12_2,chainOU.bio12_3)

combine.bio12 <- combine.chains(list_chains.bio12, burnin.prop = 0.3 )

summary.bio12 <- summary(combine.bio12)

plotSimmap.mcmc(combine.bio12, burnin = 0.3, pp.cutoff = 0.3, cex = 0.5)
plotBranchHeatMap(tree, combine.bio12, "theta", burnin = 0.3, pal = rainbow, cex = 0.5)
phenogram.density(tree, dat.bio12, burnin = 0.3, combine.bio12, pp.cutoff = 0.3)
par(mfrow=c(3, 5) )
plot(combine.bio12, auto.layout=FALSE)


