library(bayou)

tree <- read.tree("tree.tree")
mydata <- read.csv("mydata.csv")
mydata2 <- read.csv("mydata2.csv")
mydata[,1]  
setdiff(tree$tip.label, mydata[,1])
setdiff(mydata[,1],tree$tip.label)

#bio1
dat.bio1 <- mydata[,2]
names(dat.bio1) <- mydata[,1]
se.bio1 <- sqrt(mydata2[,2])
names(se.bio1) <- mydata2[,1]

priorOU.bio1 <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dk="cdpois", dtheta="dnorm"), 
                            param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                       dk=list(lambda=4, kmax=8), dsb=list(bmax=1, prob=1), 
                                       dtheta=list(mean=mean(dat.bio1), sd=1.5*sd(dat.bio1))),plot.prior = FALSE)





startpars.bio1 <- priorSim(priorOU.bio1, tree, plot=FALSE, cex=0.5)$pars[[1]]
priorOU.bio1(startpars.bio1)




set.seed(1)
#run1
mcmcOU.bio1_1 <- bayou.makeMCMC(tree, dat.bio1, SE=se.bio1, prior=priorOU.bio1, 
                               new.dir=TRUE, outname="./modelOU_rbio1_1", plot.freq=NULL) # Set up the MCMC
mcmcOU.bio1_1$run(2500000) # Run the MCMC
chainOU.bio1_1 <- mcmcOU.bio1_1$load()
summary.bio1_1 <- summary(chainOU.bio1_1)
save(chainOU.bio1_1, file = "mcmcOU.bio1_1.Rdata")

#run2
mcmcOU.bio1_2 <- bayou.makeMCMC(tree, dat.bio1, SE=se.bio1, prior=priorOU.bio1, 
                               new.dir=TRUE, outname="./modelOU_rbio1_2", plot.freq=NULL) # Set up the MCMC
mcmcOU.bio1_2$run(2500000) # Run the MCMC
chainOU.bio1_2 <- mcmcOU.bio1_2$load()
summary.bio1_2 <- summary(chainOU.bio1_2)
save(chainOU.bio1_2, file = "mcmcOU.bio1_2.Rdata")

#run3
mcmcOU.bio1_3 <- bayou.makeMCMC(tree, dat.bio1, SE=se.bio1, prior=priorOU.bio1, 
                               new.dir=TRUE, outname="./modelOU_rbio1_3", plot.freq=NULL) # Set up the MCMC
mcmcOU.bio1_3$run(2500000) # Run the MCMC
chainOU.bio1_3 <- mcmcOU.bio1_3$load()
summary.bio1_3 <- summary(chainOU.bio1_3)
save(chainOU.bio1_3, file = "mcmcOU.bio1_3.Rdata")


list_chains.bio1 <- list(chainOU.bio1_1,chainOU.bio1_2,chainOU.bio1_3)

combine.bio1 <- combine.chains(list_chains.bio1, burnin.prop = 0.3 )

summary.bio1 <- summary(combine.bio1)

plotSimmap.mcmc(combine.bio1, burnin = 0.3, pp.cutoff = 0.3, cex = 0.5)
plotBranchHeatMap(tree, combine.bio1, "theta", burnin = 0.3, pal = rainbow, cex = 0.5)
phenogram.density(tree, dat.bio1, burnin = 0.3, combine.bio1, pp.cutoff = 0.3)
par(mfrow=c(3, 5) )
plot(combine.bio1, auto.layout=FALSE)


