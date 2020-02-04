library(bayou)

tree <- read.tree("tree.tree")
mydata <- read.csv("mydata.csv")
mydata2 <- read.csv("mydata2.csv")
mydata[,1]  
setdiff(tree$tip.label, mydata[,1])
setdiff(mydata[,1],tree$tip.label)

#bio18
dat.bio18 <- mydata[,19]
names(dat.bio18) <- mydata[,1]
se.bio18 <- sqrt(mydata2[,19])
names(se.bio18) <- mydata2[,1]

priorOU.bio18 <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dk="cdpois", dtheta="dnorm"), 
                            param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                       dk=list(lambda=4, kmax=8), dsb=list(bmax=1, prob=1), 
                                       dtheta=list(mean=mean(dat.bio18), sd=1.5*sd(dat.bio18))),plot.prior = FALSE)





startpars.bio18 <- priorSim(priorOU.bio18, tree, plot=FALSE, cex=0.5)$pars[[1]]
priorOU.bio18(startpars.bio18)




set.seed(1)
#run1
mcmcOU.bio18_1 <- bayou.makeMCMC(tree, dat.bio18, SE=se.bio18, prior=priorOU.bio18, 
                               new.dir=TRUE, outname="./modelOU_rbio18_1", plot.freq=NULL) # Set up the MCMC
mcmcOU.bio18_1$run(2500000) # Run the MCMC
chainOU.bio18_1 <- mcmcOU.bio18_1$load()
summary.bio18_1 <- summary(chainOU.bio18_1)
save(chainOU.bio18_1, file = "mcmcOU.bio18_1.Rdata")

#run2
mcmcOU.bio18_2 <- bayou.makeMCMC(tree, dat.bio18, SE=se.bio18, prior=priorOU.bio18, 
                               new.dir=TRUE, outname="./modelOU_rbio18_2", plot.freq=NULL) # Set up the MCMC
mcmcOU.bio18_2$run(2500000) # Run the MCMC
chainOU.bio18_2 <- mcmcOU.bio18_2$load()
summary.bio18_2 <- summary(chainOU.bio18_2)
save(chainOU.bio18_2, file = "mcmcOU.bio18_2.Rdata")

#run3
mcmcOU.bio18_3 <- bayou.makeMCMC(tree, dat.bio18, SE=se.bio18, prior=priorOU.bio18, 
                               new.dir=TRUE, outname="./modelOU_rbio18_3", plot.freq=NULL) # Set up the MCMC
mcmcOU.bio18_3$run(2500000) # Run the MCMC
chainOU.bio18_3 <- mcmcOU.bio18_3$load()
summary.bio18_3 <- summary(chainOU.bio18_3)
save(chainOU.bio18_3, file = "mcmcOU.bio18_3.Rdata")


list_chains.bio18 <- list(chainOU.bio18_1,chainOU.bio18_2,chainOU.bio18_3)

combine.bio18 <- combine.chains(list_chains.bio18, burnin.prop = 0.3 )

summary.bio18 <- summary(combine.bio18)

plotSimmap.mcmc(combine.bio18, burnin = 0.3, pp.cutoff = 0.3, cex = 0.5)
plotBranchHeatMap(tree, combine.bio18, "theta", burnin = 0.3, pal = rainbow, cex = 0.5)
phenogram.density(tree, dat.bio18, burnin = 0.3, combine.bio18, pp.cutoff = 0.3)
par(mfrow=c(3, 5) )
plot(combine.bio18, auto.layout=FALSE)


