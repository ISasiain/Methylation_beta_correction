##set/create own home directory below:
setwd("~/Methylation/adjustBetas")

#Store the name of the individual CpG to plot
#cpg <- "cg00002719"
#cpg <- "cg09248054"
#cpg <- "cg25340711"
#cpg <- "cg16601494"
cpg <- "cg01258522"
#cpg <- "cg10244364"
#cpg <- "cg25370753"

#Loading the adjust betas environment
load("workspace_tcgaBrca_top5000.RData")

#Open files with the beta values and store them as a variable.
original_beta <- read.csv("output_betas.original.csv",sep=",",row.names=1)
corr_tum_beta <- read.csv("output_betas.tumor.csv",sep=",",row.names=1)
corr_env_beta <- read.csv("output_betas.microenvironment.csv",sep=",",row.names=1)

#Open the file with the populations detected per CpG and store it as a variable
cluster_pops <- read.csv("output_cpg.populations.csv",sep=",",row.names=1)

#Open file with slopes and intercepts
slopes_pops <- read.csv("output_reg.slopes.csv", sep=",",row.names = 1)
intercepts_pops <- read.csv("output_reg.intercepts.csv", sep=",", row.names = 1)

par(mfrow=c(1,3))

#Plot original_beta VS purity
plot(y=as.numeric(original_beta[cpg,]),
    x=1-as.numeric(purityVector), 
    col=as.numeric(cluster_pops[cpg,]),
    xlim=c(0,1), ylim=c(0,1),
    xlab="1-Purity", ylab="Beta value",
    main="Original data",
    pch=20)

#Add regression line per population
lapply(unique(as.numeric(cluster_pops[cpg,])), function(gr) 
  {
  abline(b=slopes_pops[cpg,gr], a=intercepts_pops[cpg,gr], col=gr, lwd=5)
  }
)

#Plot corrected tumor beta VS purity
#Plot original_beta VS purity
plot(y=as.numeric(corr_tum_beta[cpg,]),
     x=1-as.numeric(purityVector), 
     col=as.numeric(cluster_pops[cpg,]),
     xlim=c(0,1), ylim=c(0,1),
     xlab="Purity", ylab="Beta value",
     main="Corrected tumor betas",
     pch=20)

#Add regression line per population
lapply(unique(as.numeric(cluster_pops[cpg,])), function(gr) 
  {
  abline(h=intercepts_pops[cpg,gr], col=gr, lwd=5)
  }
)

#Plot corrected microenvironment beta VS purity.
#Plot original_beta VS purity
plot(y=as.numeric(corr_env_beta[cpg,]),
     x=1-as.numeric(purityVector), 
     col=as.numeric(cluster_pops[cpg,]),
     xlim=c(0,1), ylim=c(0,1),
     xlab="1-Purity", ylab="Beta value",
     main="Corrected microenvironment betas",
     pch=20)

#Add regression line per population
#Add regression line per population
lapply(unique(as.numeric(cluster_pops[cpg,])), function(gr) 
  {
  #The slope is added to the intercept to get the regression line of the microenvironemnt
  #The obtained value would be the same of the intercept of the beta VS purity.
  abline(h=(intercepts_pops[cpg,gr]+slopes_pops[cpg,gr]), col=gr, lwd=5)
  }
)

y <- original_beta["cg01258522",][cluster_pops["cg01258522",]==1]
x <- 1 - purityVector[cluster_pops["cg01258522",]==1]
a <- lm (y ~ x)

RSE <- summary(a)$sigma
t_val <- qt(0.25/2, 108)

to_predict <- data.frame(x=0.5)
predict(a, newdata = to_predict, interval = "prediction", level = 0.75)

predict(a, newdata = to_predict) + t_val*RSE
predict(a, newdata = to_predict) - t_val*RSE
