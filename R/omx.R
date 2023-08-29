omx <- function(data,rvEst=rep(1,ncol(data)),mubarEst=mean(data[,1]),interval=0.95,reps=500,bs.q=c(0.025,0.975),bs=TRUE) {
# Construct measurement error model and estimate parameters via FIML using OpenMx
  
dataRaw <- OpenMx::mxData( observed=data, type="raw" )
  
# residual variances
resVars <- OpenMx::mxPath( from=names(data), arrows=2,
                   free=TRUE, values=rvEst, lbound=0,
                   labels=paste("ve",1:ncol(data),sep="") )
# latent variance
latVar <- OpenMx::mxPath( from="F1", arrows=2, lbound=0,
                  free=TRUE, values=var(data[,1]), labels ="sigma2" )
# factor loadings
facLoads <- OpenMx::mxPath( from="F1", to=names(data), arrows=1,
                    free=c(FALSE,rep(TRUE,ncol(data)-1)), values=rep(1,ncol(data)),
                    labels = paste("b",1:ncol(data),sep=""))
# means
means <- OpenMx::mxPath( from="one", to=c(names(data),"F1"),
                 free=c(FALSE,rep(TRUE,ncol(data))), values=c(rep(0,ncol(data)),mubarEst),
                 labels =c(paste("a",1:ncol(data),sep=""),"mubar") )

# Create algebras to take square roots of variances for convenience
alg1 <- OpenMx::mxAlgebra(sqrt(sigma2),"sigma")
alg2 <- OpenMx::mxAlgebra(sqrt(ve1),"se1")
ALGse <- paste("algse",2:ncol(data)," <- OpenMx::mxAlgebra(sqrt(ve",2:ncol(data),")/b",2:ncol(data),",'base",2:ncol(data),"')",sep="")

for(i in 1:ncol(data)){ eval(parse(text=ALGse[i])) }

lst <- c(paste("b",2:ncol(data),sep=""),"se1",paste("base",2:ncol(data),sep=""),paste("a",2:ncol(data),sep=""),"mubar","sigma")
ci <- mxCI(lst,interval=interval)

# Put model together from individual pieces shown above
run <- paste("oneFactorModel <- OpenMx::mxModel('Common Factor Model Path Specification', type='RAM',manifestVars=c(names(data)), latentVars='F1',dataRaw, resVars, latVar, facLoads, means, alg1, alg2,", paste("algse",2:ncol(data),collapse=",",sep=""),", ci)")
eval(parse(text=run))

# Fit the model and estimate parameters
oneFactorFit <- OpenMx::mxTryHard(oneFactorModel,intervals=TRUE)

CI <- summary(oneFactorFit)$CI
rownames(CI) <- lst

list(fit=oneFactorFit,ci=CI,abs=alpha.beta.sigma(summary(oneFactorFit)$parameters[,c(1,5,6)]),model=oneFactorModel)

# Boostrap by default
if(bs) {
vboot <- OpenMx::mxBootstrap(oneFactorFit,replications=reps)

# Create a vector of character strings of R code to be executed 
rw <- paste("r",1:(3*ncol(data))," <- OpenMx::mxBootstrapEval(",lst,",vboot,bq=c(",paste(bs.q,collapse=",",sep=""),"))",sep="") 
  
# Execute each element of vector
for(i in 1:(3*ncol(data))) {eval(parse(text=rw[i]))}
  
# Combine and name all the results for return
boot.q <- eval(parse(text=paste("rbind(",paste("r",1:(3*ncol(data)),collapse=",",sep=""),")",sep="")))
rownames(boot.q) <- lst
  
# Return the results
list(fit=oneFactorFit,ci=CI,boot=vboot,q.boot=boot.q,bsReps=reps,abs=alpha.beta.sigma(summary(oneFactorFit)$parameters[,c(1,5,6)]),model=oneFactorModel)
}
  
}