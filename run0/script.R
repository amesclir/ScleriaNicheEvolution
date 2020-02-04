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
mcmcOU.bio1 <- bayou.makeMCMC(tree, dat.bio1, SE=se.bio1, prior=priorOU.bio1, 
                               new.dir=TRUE, outname="./modelOU_rbio1", plot.freq=NULL) # Set up the MCMC
mcmcOU.bio1$run(2500000) # Run the MCMC
chainOU.bio1 <- mcmcOU.bio1$load()

save(chainOU.bio1, file = "mcmcOU.bio1.Rdata")

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
mcmcOU.bio4 <- bayou.makeMCMC(tree, dat.bio4, SE=se.bio4, prior=priorOU.bio4, 
                               new.dir=TRUE, outname="./modelOU_rbio4", plot.freq=NULL) # Set up the MCMC
mcmcOU.bio4$run(2500000) # Run the MCMC
chainOU.bio4 <- mcmcOU.bio4$load()

save(chainOU.bio4, file = "mcmcOU.bio4.Rdata")


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
mcmcOU.bio7 <- bayou.makeMCMC(tree, dat.bio7, SE=se.bio7, prior=priorOU.bio7, 
                               new.dir=TRUE, outname="./modelOU_rbio7", plot.freq=NULL) # Set up the MCMC
mcmcOU.bio7$run(2500000) # Run the MCMC
chainOU.bio7 <- mcmcOU.bio7$load()

save(chainOU.bio7, file = "mcmcOU.bio7.Rdata")

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
mcmcOU.bio12 <- bayou.makeMCMC(tree, dat.bio12, SE=se.bio12, prior=priorOU.bio12, 
                               new.dir=TRUE, outname="./modelOU_rbio12", plot.freq=NULL) # Set up the MCMC
mcmcOU.bio12$run(2500000) # Run the MCMC
chainOU.bio12 <- mcmcOU.bio12$load()

save(chainOU.bio12, file = "mcmcOU.bio12.Rdata")

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
mcmcOU.bio15 <- bayou.makeMCMC(tree, dat.bio15, SE=se.bio15, prior=priorOU.bio15, 
                               new.dir=TRUE, outname="./modelOU_rbio15", plot.freq=NULL) # Set up the MCMC
mcmcOU.bio15$run(2500000) # Run the MCMC
chainOU.bio15 <- mcmcOU.bio15$load()

save(chainOU.bio15, file = "mcmcOU.bio15.Rdata")


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
mcmcOU.bio18 <- bayou.makeMCMC(tree, dat.bio18, SE=se.bio18, prior=priorOU.bio18, 
                               new.dir=TRUE, outname="./modelOU_rbio18", plot.freq=NULL) # Set up the MCMC
mcmcOU.bio18$run(2500000) # Run the MCMC
chainOU.bio18 <- mcmcOU.bio18$load()

save(chainOU.bio18, file = "mcmcOU.bio18.Rdata")
