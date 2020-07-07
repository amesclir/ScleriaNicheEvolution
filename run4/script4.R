library(bayou)

tree <- read.tree("tree.tree")
mydata <- read.csv("mydata.csv")
mydata2 <- read.csv("mydata2.csv")
mydata[,1]  
setdiff(tree$tip.label, mydata[,1])
setdiff(mydata[,1],tree$tip.label)




#bio4
dat.bio4 <- mydata[,5]
names(dat.bio4) <- mydata[,1]
se.bio4 <- sqrt(mydata2[,5])
names(se.bio4) <- mydata2[,1]

priorOU.bio4 <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dk="cdpois", dtheta="dnorm"), 
                            param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                       dk=list(lambda=4, kmax=8), dsb=list(bmax=1, prob=1), 
                                       dtheta=list(mean=mean(dat.bio4), sd=1.5*sd(dat.bio4))),plot.prior = FALSE)





startpars.bio4 <- priorSim(priorOU.bio4, tree, plot=FALSE, cex=0.5)$pars[[1]]
priorOU.bio4(startpars.bio4)




set.seed(1)
#run1
mcmcOU.bio4_1 <- bayou.makeMCMC(tree, dat.bio4, SE=se.bio4, prior=priorOU.bio4, 
                               new.dir=TRUE, outname="./modelOU_rbio4_1", plot.freq=NULL) # Set up the MCMC
mcmcOU.bio4_1$run(2500000) # Run the MCMC
chainOU.bio4_1 <- mcmcOU.bio4_1$load()
summary.bio4_1 <- summary(chainOU.bio4_1)
save(chainOU.bio4_1, file = "mcmcOU.bio4_1.Rdata")

#run2
mcmcOU.bio4_2 <- bayou.makeMCMC(tree, dat.bio4, SE=se.bio4, prior=priorOU.bio4, 
                               new.dir=TRUE, outname="./modelOU_rbio4_2", plot.freq=NULL) # Set up the MCMC
mcmcOU.bio4_2$run(2500000) # Run the MCMC
chainOU.bio4_2 <- mcmcOU.bio4_2$load()
summary.bio4_2 <- summary(chainOU.bio4_2)
save(chainOU.bio4_2, file = "mcmcOU.bio4_2.Rdata")

#run3
mcmcOU.bio4_3 <- bayou.makeMCMC(tree, dat.bio4, SE=se.bio4, prior=priorOU.bio4, 
                               new.dir=TRUE, outname="./modelOU_rbio4_3", plot.freq=NULL) # Set up the MCMC
mcmcOU.bio4_3$run(2500000) # Run the MCMC
chainOU.bio4_3 <- mcmcOU.bio4_3$load()
summary.bio4_3 <- summary(chainOU.bio4_3)
save(chainOU.bio4_3, file = "mcmcOU.bio4_3.Rdata")


list_chains.bio4 <- list(chainOU.bio4_1,chainOU.bio4_2,chainOU.bio4_3)

combine.bio4 <- combine.chains(list_chains.bio4, burnin.prop = 0.3 )

summary.bio4 <- summary(combine.bio4)

plotSimmap.mcmc(combine.bio4, burnin = 0.3, pp.cutoff = 0.3, cex = 0.5)
plotBranchHeatMap(tree, combine.bio4, "theta", burnin = 0.3, pal = rainbow, cex = 0.5)
phenogram.density(tree, dat.bio4, burnin = 0.3, combine.bio4, pp.cutoff = 0.3)
par(mfrow=c(3, 5) )
plot(combine.bio4, auto.layout=FALSE)


