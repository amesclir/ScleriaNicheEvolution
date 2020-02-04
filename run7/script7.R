library(bayou)

tree <- read.tree("tree.tree")
mydata <- read.csv("mydata.csv")
mydata2 <- read.csv("mydata2.csv")
mydata[,1]  
setdiff(tree$tip.label, mydata[,1])
setdiff(mydata[,1],tree$tip.label)




#bio7
dat.bio7 <- mydata[,8]
names(dat.bio7) <- mydata[,1]
se.bio7 <- sqrt(mydata2[,8])
names(se.bio7) <- mydata2[,1]

priorOU.bio7 <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dk="cdpois", dtheta="dnorm"), 
                            param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1),
                                       dk=list(lambda=4, kmax=8), dsb=list(bmax=1, prob=1), 
                                       dtheta=list(mean=mean(dat.bio7), sd=1.5*sd(dat.bio7))),plot.prior = FALSE)





startpars.bio7 <- priorSim(priorOU.bio7, tree, plot=FALSE, cex=0.5)$pars[[1]]
priorOU.bio7(startpars.bio7)




set.seed(1)
#run1
mcmcOU.bio7_1 <- bayou.makeMCMC(tree, dat.bio7, SE=se.bio7, prior=priorOU.bio7, 
                               new.dir=TRUE, outname="./modelOU_rbio7_1", plot.freq=NULL) # Set up the MCMC
mcmcOU.bio7_1$run(2500000) # Run the MCMC
chainOU.bio7_1 <- mcmcOU.bio7_1$load()
summary.bio7_1 <- summary(chainOU.bio7_1)
save(chainOU.bio7_1, file = "mcmcOU.bio7_1.Rdata")

#run2
mcmcOU.bio7_2 <- bayou.makeMCMC(tree, dat.bio7, SE=se.bio7, prior=priorOU.bio7, 
                               new.dir=TRUE, outname="./modelOU_rbio7_2", plot.freq=NULL) # Set up the MCMC
mcmcOU.bio7_2$run(2500000) # Run the MCMC
chainOU.bio7_2 <- mcmcOU.bio7_2$load()
summary.bio7_2 <- summary(chainOU.bio7_2)
save(chainOU.bio7_2, file = "mcmcOU.bio7_2.Rdata")

#run3
mcmcOU.bio7_3 <- bayou.makeMCMC(tree, dat.bio7, SE=se.bio7, prior=priorOU.bio7, 
                               new.dir=TRUE, outname="./modelOU_rbio7_3", plot.freq=NULL) # Set up the MCMC
mcmcOU.bio7_3$run(2500000) # Run the MCMC
chainOU.bio7_3 <- mcmcOU.bio7_3$load()
summary.bio7_3 <- summary(chainOU.bio7_3)
save(chainOU.bio7_3, file = "mcmcOU.bio7_3.Rdata")


list_chains.bio7 <- list(chainOU.bio7_1,chainOU.bio7_2,chainOU.bio7_3)

combine.bio7 <- combine.chains(list_chains.bio7, burnin.prop = 0.3 )

summary.bio7 <- summary(combine.bio7)

plotSimmap.mcmc(combine.bio7, burnin = 0.3, pp.cutoff = 0.3, cex = 0.5)
plotBranchHeatMap(tree, combine.bio7, "theta", burnin = 0.3, pal = rainbow, cex = 0.5)
phenogram.density(tree, dat.bio7, burnin = 0.3, combine.bio7, pp.cutoff = 0.3)
par(mfrow=c(3, 5) )
plot(combine.bio7, auto.layout=FALSE)


