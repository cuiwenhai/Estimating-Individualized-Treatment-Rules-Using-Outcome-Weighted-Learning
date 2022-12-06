


library(quadprog)
library(Dykstra)




GetDesign <- function(n, p=50, 
                      unif.min=-1, unif.max=1,
                      debugging.DI=FALSE) {
  # generate (X, Z), the patient characteristics that aren't associated with treatment or prognosis
  X <- matrix(0, nrow=n, ncol=p)
  for (j in 1:p) {
    X[, j] <- runif(n=n, unif.min, unif.max)
  }
  #colnames(X) <- paste0("input_", 1:p)
  T <- as.matrix(rbinom(n=n, size=1, prob=0.5))
  T<-2*T-1
  return(list("X"=X, "T"=T))
}


Expit <- function(x) exp(x) / (1 + exp(x))

GetResponse <- function(scenario,
                        n, T, X,
                        shift.for.continuous=c("none", "min.and.unit.sd", "min", "10", "50"), 
                        epsilon.mu=0, epsilon.sigma=1) {
  my.shift.for.continuous <- match.arg(shift.for.continuous)
  #stopifnot(n == nrow(X))
  # compute ``score'' underlying response
  S <- rep(NA, n)
  S.mat <- matrix(nrow=n, ncol=3)
  L=rep(0,n)
  if (scenario %in% c("default", "default.continuous", "noise.added", "noise.added.continuous")) {
    S[T == -1] <- 1+2*X[T == -1,1]+X[T == -1,2]+0.5*X[T == -1,3]+(1-X[T == -1,1]-X[T == -1,2])*(-1)
    S[T == 1] <- 1+2*X[T==1,1]+X[T==1,2]+0.5*X[T==1,3]+(1-X[T==1,1]-X[T==1,2])*1
  } else if (scenario %in% "chapter3") {
    S[T == 0] <- 0 + X[T==0, 1] * 0 + L[T==0, 1] * 1.5 + X[T==0, 1] * -0.55 + 0 * 0 + 0 * 0
    S[T == 1] <- -1.4 + X[T==1, 1] * 0 + L[T==1, 1] * 1.5 + X[T==1, 1] * 0.55 + 0 * 0 + 0 * 0
  } else if (scenario %in% c("main.effect.added", "main.effect.added.continuous")) {
    S[T == 0] <- 0 + X[T==0, 1] * 1 + L[T==0, 1] * 1.5 + X[T==0, 1] * -0.7 + 0 * 0 + 0 * 0
    S[T == 1] <- -1.4 + X[T==1, 1] * 1 + L[T==1, 1] * 1.5 + X[T==1, 1] * 0.4 + 0 * 0 + 0 * 0
  } else if (scenario %in% c("quadratic.interaction.added", "quadratic.interaction.added.continuous")) {
    S[T == 0] <- -0.2 + X[T==0, 1] * 0 + L[T==0, 1] * 1.5 + X[T==0, 1] * 1 + (X[T==0, 1]^2) * -1 + 0 * 0
    S[T == 1] <- -0.5 + X[T==1, 1] * 0 + L[T==1, 1] * 1.5 + X[T==1, 1] * -0.3 + (X[T==1, 1]^2) * 0.8 + 0 * 0
  } else if (scenario %in% c("piecewise.cubic.interaction.added", "piecewise.cubic.interaction.added.continuous")) {
    S[T == 0] <- 0 + X[T==0, 1] * 0 + L[T==0, 1] * 1.5 + X[T==0, 1] * -0.7 + (X[T==0, 1]^2) * -0.3 + 0
    S[T == 1] <- -1 + X[T==1, 1] * 0 + L[T==1, 1] * 1.5 + (X[T==1, 1] < 0.8) * (X[T==1, 1]  * 0.4  + X[T==1, 1]^3 * 2) +
      (X[T==1, 1] > 0.8) * (X[T==1, 1]  * 2.3 + X[T==1, 1]^3 * -1)
  } else if (scenario %in% c("high.dimensional.noise.added", "high.dimensional.noise.added.continuous")) {
    S[T == 0] <- 0 + X[T==0, 1] * 0 + L[T==0, 1] * 1.5 + X[T==0, 1] * -0.7 + 0 * 0 + 0 * 0
    S[T == 1] <- -1.4 + X[T==1, 1] * 0 + L[T==1, 1] * 1.5 + X[T==1, 1] * 0.4 + 0 * 0 + 0 * 0
  } else if (scenario %in% c("high.dimensional.noise.and.prognostic.added", "high.dimensional.noise.and.prognostic.added.continuous")) {
    S[T == 0] <- 0 + X[T==0, 1] * 0 + L[T==0, 1] * 1.5 + X[T==0, 1] * -0.7 + 0 * 0 + 0 * 0 +
      X[T==0, 2] * 0.2 + X[T==0, 3] * 0.3 + X[T==0, 4] * -0.3 + X[T==0, 5] * -0.2 + X[T==0, 6] * 0.4
    S[T == 1] <- -1.4 + X[T==1, 1] * 0 + L[T==1, 1] * 1.5 + X[T==1, 1] * 0.4 + 0 * 0 + 0 * 0 +
      X[T==1, 2] * 0.2 + X[T==1, 3] * 0.3 + X[T==1, 4] * -0.3 + X[T==1, 5] * -0.2 + X[T==1, 6] * 0.4
  } else if (scenario %in% c("default.piecewise.constant", "default.piecewise.constant.continuous",
                             "two.way.interactions.piecewise.constant", "two.way.interactions.piecewise.constant.continuous",
                             "three.way.interactions.piecewise.constant", "three.way.interactions.piecewise.constant.continuous",
                             "nonlinear.interactions.piecewise.constant", "nonlinear.interactions.piecewise.constant.continuous")) {
    S.mat[T==0, 1] <- 0.8 * (X[T==0, 1] > 0 & X[T==0, 1] < 0.7) + -1 * (X[T==0, 1] > 0.7 & X[T==0, 1] < 1.5) + 1 * (X[T==0, 1] > 1.5 & X[T==0, 1] < 2)
    S.mat[T==0, 2] <- 0.3 * (X[T==0, 2] > 0 & X[T==0, 2] < 0.5) + -0.7 * (X[T==0, 2] > 0.5 & X[T==0, 2] < 1) + 0 * (X[T==0, 2] > 1 & X[T==0, 2] < 1.3) + -0.5 * (X[T==0, 2] > 1.3 & X[T==0, 2] < 2)
    S.mat[T==0, 3] <- -0.3 * (X[T==0, 3] > 0 & X[T==0, 3] < 0.8) + 0.4 * (X[T==0, 3] > 0.8 & X[T==0, 3] < 1.3) + -0 * (X[T==0, 3] > 1.3 & X[T==0, 3] < 2)
    
    S.mat[T==1, 1] <- 0.5 * (X[T==1, 1] > 0 & X[T==1, 1] < 0.5) + 0 * (X[T==1, 1] > 0.5 & X[T==1, 1] < 1.3) + -0.5 * (X[T==1, 1] > 1.3 & X[T==1, 1] < 2)
    S.mat[T==1, 2] <- -0.3 * (X[T==1, 2] > 0 & X[T==1, 2] < 0.6) + -1 * (X[T==1, 2] > 0.6 & X[T==1, 2] < 1.1) + 0.6 * (X[T==1, 2] > 1.1 & X[T==1, 2] < 1.4) + 0.2 * (X[T==1, 2] > 1.4 & X[T==1, 2] < 2)
    S.mat[T==1, 3] <- 0.5 * (X[T==1, 3] > 0 & X[T==1, 3] < 0.4) + -1 * (X[T==1, 3] > 0.4 & X[T==1, 3] < 0.9) + 1.1 * (X[T==1, 3] > 0.9 & X[T==1, 3] < 2)
    if (scenario %in% c("default.piecewise.constant", "default.piecewise.constant.continuous")) {
      S[T==0] <- apply(S.mat[T==0, ], 1, sum) + 0 + L[T==0, 1] * 0.5
      S[T==1] <- apply(S.mat[T==1, ], 1, sum) + 0 + L[T==1, 1] * 0.5
    } else if (scenario %in% c("two.way.interactions.piecewise.constant", "two.way.interactions.piecewise.constant.continuous")) {
      S[T==0] <- S.mat[T==0, 1] + S.mat[T==0, 2] + S.mat[T==0, 3] +
        -0.2 * (S.mat[T==0, 1] * S.mat[T==0, 2]) +
        0.2 * (S.mat[T==0, 1] * S.mat[T==0, 3]) +
        -0.3 * (S.mat[T==0, 2] * S.mat[T==0, 3]) +
        0 + L[T==0, 1] * 0.5
      S[T==1] <- S.mat[T==1, 1] + S.mat[T==1, 2] + S.mat[T==1, 3] +
        0.1 * (S.mat[T==1, 1] * S.mat[T==1, 2]) +
        0.3 * (S.mat[T==1, 1] * S.mat[T==1, 3]) +
        -0.1 * (S.mat[T==1, 2] * S.mat[T==1, 3]) +
        0 + L[T==1, 1] * 0.5
    } else if (scenario %in% c("three.way.interactions.piecewise.constant", "three.way.interactions.piecewise.constant.continuous")) {
      S[T==0] <- S.mat[T==0, 1] + S.mat[T==0, 2] + S.mat[T==0, 3] +
        -0.5 * (S.mat[T==0, 1] * S.mat[T==0, 2] * S.mat[T==0, 3]) +
        0 + L[T==0, 1] * 0.5
      S[T==1] <- S.mat[T==1, 1] + S.mat[T==1, 2] + S.mat[T==1, 3] +
        0.5 * (S.mat[T==1, 1] * S.mat[T==1, 2] * S.mat[T==1, 3]) +
        0 + L[T==1, 1] * 0.5
    } else if (scenario %in% c("nonlinear.interactions.piecewise.constant", "nonlinear.interactions.piecewise.constant.continuous")) {
      S[T==0] <- S.mat[T==0, 1] + S.mat[T==0, 2] + S.mat[T==0, 3] +
        -0.3 * (S.mat[T==0, 1] > 1) * (S.mat[T==0, 2]) ^ 2 +
        0.4 * (S.mat[T==0, 2] > 0.5) * (S.mat[T==0, 3]) ^ 2 +
        -0.1 * (S.mat[T==0, 1] < 1.5) * S.mat[T==0, 3] ^3
      0 + L[T==0, 1] * 0.5
      S[T==1] <- S.mat[T==1, 1] + S.mat[T==1, 2] + S.mat[T==1, 3] +
        0.4 * (S.mat[T==1, 1] > 1) * (S.mat[T==1, 2]) ^ 2 +
        -0.5 * (S.mat[T==1, 2] > 0.5) * (S.mat[T==1, 3]) ^ 2 +
        -0.6 * (S.mat[T==1, 1] < 1.5) * S.mat[T==1, 3] ^3
      0 + L[T==0, 1] * 0.5
    }
  }
  # convert score into probability of response (for binary outcome)
  prob.binary.response <- Expit(S)
  # get response variables/..
  Y.binary <- rbinom(n, size=1, prob=prob.binary.response)
  if (my.shift.for.continuous %in% c(10, 50)) {
    Y.continuous <- S + rnorm(n, mean=epsilon.mu, sd=epsilon.sigma) + as.numeric(my.shift.for.continuous)
  } else if (my.shift.for.continuous %in% c("none")) {
    Y.continuous <- S + rnorm(n, mean=epsilon.mu, sd=epsilon.sigma)
  } else if (my.shift.for.continuous %in% c("min")) {
    Y.continuous.first <- S + rnorm(n, mean=epsilon.mu, sd=epsilon.sigma)
    Y.continuous <- Y.continuous.first + abs(min(Y.continuous.first)) * 1.01
  } else if (my.shift.for.continuous %in% c("min.and.unit.sd")) {
    Y.continuous.first <- S + rnorm(n, mean=epsilon.mu, sd=epsilon.sigma)
    Y.continuous <- Y.continuous.first + abs(min(Y.continuous.first)) * 1.01
    Y.continuous <- Y.continuous / sd(Y.continuous)
  }
  return(list("S"=S, "S.mat"=S.mat,
              "prob.binary.response"=prob.binary.response,
              "Y.binary"=Y.binary, "Y.continuous"=Y.continuous))
}




#################simulate data########################
n=100
X=as.matrix(GetDesign(n)$X)
T=GetDesign(n)$T
Y=GetResponse(scenario="default",n,T,X,shift.for.continuous=c("none"), 
              epsilon.mu=0, epsilon.sigma=1)

Y=as.matrix(Y$Y.continuous,n,1)

p=ncol(X)
#################求解model参数########################
sgn<-function(x)
{ if(x>0)
{return(1)}else{return(-1)}
}###符号函数

decision.function<-function(x,X,T,Alpha.solution,bete0.hat,n){
  k.new=matrix(0,n,1)
 
  for (i in 1:n) {
    
    a=matrix(X[i,],50,1)
    b=matrix(x,50,1)
    k.new[i,]=kenel(a,b)
  }
  k<-k.new
  d<-t(Alpha.solution*T)%*%k.new+bete0.hat
  return(d) 
}##决策函数


kenel<-function(X,Y){ 
  sigma=1
  k=exp(-(t(X-Y)%*%(X-Y))*(sigma^2))
return(k)}###核函数
#kenel<-function(X,Y){ 
#sigma=1
#k=t(X)%*%(Y)
#return(k)}


model<-function(X,T,Y){###模型求解参数
n=nrow(X)
p=ncol(X)
k=matrix(0,n,n)
lambda<-1
kapa<-0.5/lambda
e<-rep(1,n)
###核矩阵
#for (i in 1:n) {
#  for (j in 1:n) { 
#    a=matrix(X[i,],p,1)
#    b=matrix(X[j,],p,1)
#    k[i,j]=kenel(a,b)*T[i,]*T[j,]
#  }}
for (i in 1:n) {
  for (j in 1:n) { 
    a=matrix(X[i,],p,1)
    b=matrix(X[j,],p,1)
    k[i,j]=kenel(a,b)
  }}
k<-k++0.01*diag(e)
###定义损失函数
Alpha<-matrix(0,n,1)#要求解的参数
d<-(matrix(rep(2,n),n,1))/T
loss<-function(T,X,Y,Alpha){
  loss=(t(Alpha)%*%k%*%(Alpha)-t(Alpha)%*%(d))*0.5
  return(loss)
  }
###根据损失函数求解二次规划问题Alpha

Control.matrix<--rbind(rep(1,n),-diag(as.vector(T),n,n),diag(as.vector(T),n,n))
#Control.matrix<--rbind(t(T),-diag(e),diag(e))
bound.vec.1<-c(0,(Y/0.5)*kapa,rep(0,n))
for (i in 1:n+1) {
 if( bound.vec.1[i]<0){
   bound.vec.1[i]=0.001
    }
  }
bound.vec<--bound.vec.1###控制条件
### 求解参数
Alpha<-dykstra(Dmat=k, dvec=d, Amat=t(Control.matrix), bvec=bound.vec, meq = 1, factorized = FALSE,
               maxit = NULL, eps = NULL)
#Alpha.solution.1<-Alpha$solution
Alpha.solution.1<-Alpha$solution/T

Alpha.solution<-as.matrix(Alpha.solution.1,n,1)

length(bound.vec)
ncol(t(Control.matrix))
bound.vec<--bound.vec
#Alpha<-solve.QP(Dmat=k, dvec=d, Amat=t(Control.matrix), bvec=bound.vec, meq=1, factorized=FALSE)
#t(T)%*%Alpha.solution         ###应该约等于0
#cbind(bound.vec,rbind(0,Alpha.solution),rbind(0,T))
#loss(T,X,Y,Alpha.solution)
######求解决策函数的参数bete0

bete0.hat<-0
k1<-0
bete0.hat1<-matrix(0,n,1)
for (i in 1:n) {
 if( 0<Alpha.solution[i,]& Alpha.solution[i,]<bound.vec[i+1]-0.01){
   k1=k1+1
   bete0.hat1[i,1]<-T[i]-decision.function(x=X[i,],X,T,Alpha.solution,bete0.hat=0,n)
   bete0.hat<-bete0.hat1[i,1]/k1+bete0.hat*(1-1/k1)
 }
 
}
#cbind(rbind(0,bete0.hat1) ,rbind(0,T),rbind(0,Alpha.solution),bound.vec)  
bete0.hat

#bete.hat<-(t(X)%*%(Alpha.solution*T))
#bete.hat

return(list("Alpha.solution"=Alpha.solution,"bete0.hat"=bete0.hat,"p"=p))

}


model.0<-model(X,T,Y)
Alpha.solution.0<-model.0[[1]]
bete0.hat.0<-model.0[[2]]
p<-model.0[[3]]

bete0.hat.0

###########生成测试集

test.number<-1000
x.test<- matrix(0, nrow=test.number, ncol=p)
test.number<-1000
x.test<- matrix(0, nrow=test.number, ncol=p)
for (j in 1:test.number) {
  x.test[j, ] <- runif(n=p, -1, 1)}##产生测试集

###########原模型进行测试，查看在测试集上表现
misclass<-0##开始测试
mse<-matrix(0,10,1)
for (j in 1:test.number) {
  decision=sgn(decision.function(x.test[j, ],X,T,Alpha.solution.0,bete0.hat.0,n))
  a=x.test[j,c(1,2)]
  
  if( (sum(a)<1&decision==-1) | (sum(a)>1& decision==1) ){misclass=misclass+1
  }
print(j)
  }
mse[1,]<-misclass/test.number
mse



################################################
#####################boosting改变数据的分布.B次，发现没有明显效果
################################################



change.distribution<-function(X,T,Y,Alpha.solution.0,bete0.hat.0){
  n<-nrow(X)
  prob<-matrix(0,n,1)
  for (j in 1:n) {
    #prob[j]=decision.function(X[j,],X,T,Alpha.solution.0,bete0.hat.0,n)
    prob[j]=Y[j]
    print(j)
  }
  omega<-1
  prob.1<-exp(omega*abs(prob))
  prob.real<-prob.1/sum(prob.1)
  sample<-sample(1:n,size = n,replace=TRUE,prob=prob.real)
  plot(density(sample))
  X.1<-matrix(0,n,p)
  T.1<-matrix(0,n,1)
  Y.1<-matrix(0,n,1)
  for (j in 1:n) {
    s<-sample[j]
    X.1[j,]<-X[s,]
    T.1[j,]<-T[s,]
    Y.1[j,]<-Y[s,]
  }
  resulat<-list(X.1,T.1,Y.1)
  return(resulat)
}
B=5
mse<-matrix(0,B,1)
X.boostig<-X
T.boostig<-T
Y.boostig<-Y

Alpha.solution.B<-matrix(0,n,B)
bete0.hat.B<-matrix(0,1,B)

Alpha.solution.B[,1]<-Alpha.solution.0
bete0.hat.B[,1]<-bete0.hat.0

for (b in 1:B) {
  new.data<-change.distribution(X.boostig,T.boostig,Y.boostig,Alpha.solution.B[,b],bete0.hat.B[,b])
  X.boostig<-new.data[[1]]
  T.boostig<-new.data[[2]]
  Y.boostig<-new.data[[3]]
  ###########boosting再次进行求解
  model.1<-model(X.boostig,T.boostig,Y.boostig)
  Alpha.solution.B[,b]<-model.1[[1]]
  bete0.hat.B[,b]<-model.1[[2]]
}
###########测试集查看效果
mse<-matrix(0,B,1)
for (l in 1:5) {
  misclass<-0
  for (j in 1:test.number) {
    decision.vector<-matrix(0,B,1)
    for (b in 1:l) {
      decision.vector[b,]<-decision.function(x.test[j, ],X.boostig,T.boostig,Alpha.solution.B[,b],bete0.hat.B[,b],n)
    }
    decision=sgn( sum(b))
    a=x.test[j,c(1,2)]
    
    if( (sum(a)<1&decision==-1) | (sum(a)>1& decision==1) ){misclass=misclass+1
    }
    print(j)
  }
  mse[l,]<-misclass/test.number
}

mse

#####################误分类数据再学习，发现并没有提升效果
X.mis<-matrix(0,n,p)
T.mis<-matrix(0,n,1)
Y.mis<-matrix(0,n,1)
i=0
for (j in 1:n) {
   
  decision=sgn(decision.function(X[j, ],X,T,Alpha.solution.0,bete0.hat.0,n))
  a=X[j,c(1,2)]
  
  if( (sum(a)<1&decision==-1) | (sum(a)>1& decision==1) ){
    X.mis[i,]<-X[j,]
    T.mis[i,]<-T[j,]
    Y.mis[i,]<-Y[j,]
    i=i+1
  }
  print(j)
}
X.mis.real<-X.mis[-14:-n,]
T.mis.real<-as.matrix(T[-14:-n,],14,1)
Y.mis.real<-as.matrix(Y[-14:-n,],14,1)


model.1<-model(X.mis.real,T.mis.real,Y.mis.real)
Alpha.solution.1<-model.1[[1]]
bete0.hat.1<-model.1[[2]]
p<-model.1[[3]]
###测试

misclass<-0  ##开始测试
mse<-matrix(0,10,1)
for (j in 1:test.number) {
 
decision=sgn( decision.function(x.test[j, ],X,T,Alpha.solution.0,bete0.hat.0,n)
                            +decision.function(x.test[j, ],X.mis.real,T.mis.real,Alpha.solution.1,bete0.hat.1,13) ) 
  a=x.test[j,c(1,2)]
  
  if( (sum(a)<1&decision==-1) | (sum(a)>1& decision==1) ){misclass=misclass+1
  }
  print(j)
}
mse[1,]<-misclass/test.number
misclass<-0 
for (j in 1:test.number) {
  
  decision=sgn( decision.function(x.test[j, ],X.mis.real,T.mis.real,Alpha.solution.1,bete0.hat.1,13) ) 
  a=x.test[j,c(1,2)]
  
  if( (sum(a)<1&decision==-1) | (sum(a)>1& decision==1) ){misclass=misclass+1
  }
  print(j)
}
mse[2,]<-misclass/test.number
misclass<-0 
for (j in 1:test.number) {
  
  decision=sgn( decision.function(x.test[j, ],X,T,Alpha.solution.0,bete0.hat.0,n)) 
  a=x.test[j,c(1,2)]
  
  if( (sum(a)<1&decision==-1) | (sum(a)>1& decision==1) ){misclass=misclass+1
  }
  print(j)
}
mse[3,]<-misclass/test.number
mse













#########################  结束   #############################################

#x<-matrix(c(-1,-1,rep(0,48)),50,1)

########验证决策
Y.no<-matrix(0,100,1)
Y.yes<-matrix(0,100,1)
for (i in 1:100) {
 Y.x.treament=GetResponse(scenario="default",n=1,T=1,t(x),shift.for.continuous=c("none"), 
                                epsilon.mu=0, epsilon.sigma=1)
 Y.yes[i,1]=Y.x.treament$Y.continuous
 Y.x.no.treament=GetResponse(scenario="default",n=1,T=-1,t(x),shift.for.continuous=c("none"), 
                                   epsilon.mu=0, epsilon.sigma=1)
 Y.no[i,1]=Y.x.no.treament$Y.continuous
}
plot(c(1:100),Y.yes,col="red",type="l",ylim=c(-5,8) )
lines(c(1:100),Y.no,col="blue",type="l")  

mean(Y.yes)
mean(Y.no)
sgn(decision)
decision
x[c(1,2),]



















###########mep测试表示mep及其之前的是等式约束，mep以后是不等式约束。
## Assume we want to minimize: -(0 1 ) %*% b + 1/2 b^T b
## under the constraints:      A^T b >= b0
## with b0 = (1,1)^T
## and      (1  0  ) 
##      A = (0  1 )
##          
## we can use solve.QP as follows:
##
Dmat       <- matrix(0,2,2)
diag(Dmat) <- c(1,1)
dvec       <- c(0,1)
Amat       <- matrix(c(1,0,0,1),2,2)
bvec       <- as.matrix(c(-1,0.5),2,1)
solve.QP(Dmat,dvec,Amat,bvec=t(bvec),meq=1)
#########