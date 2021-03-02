GGRidge<-function(data, kg=5, sq=c(0.01,5,0.01), k=5, PE=50){
  lm.ridge.univariate<-function (x, y, lambda = 0, scale = TRUE) {
    r <- length(lambda)
    xs <- scale(x, scale = scale)
    ys <- scale(y, scale = scale)
    sx <- 1
    sy <- 1
    if (scale == TRUE) {
      sx <- sd(x)
      sy <- sd(y)
    }
    b <- sum(xs * ys)/(sum(xs^2) + lambda)
    b <- b * sy/sx
    inter <- mean(y) - b * mean(x)
    coefficients <- cbind(inter, b)
    return(coefficients)
  }
  ridge.cv<-function (X, y, lambda = NULL, scale = TRUE, k = 10, plot.it = FALSE) {
   
    if (is.vector(X) == TRUE) {
      X <- matrix(X, ncol = 1)
    }
    if (is.null(lambda) == TRUE) {
      ss <- seq(-10, -1, length = 1000)
      ss <- 10^ss
      n <- nrow(X)
      nn <- n - floor(n/k)
      lambda <- ss * nn * ncol(X)
    }
    cv <- rep(0, length(lambda))
    n <- nrow(X)
    all.folds <- split(sample(1:n), rep(1:k, length = n))
    for (i in seq(k)) {
      omit <- all.folds[[i]]
      Xtrain = X[-omit, , drop = FALSE]
      ytrain = y[-omit]
      Xtest = X[omit, , drop = FALSE]
      ytest = y[omit]
      if ((is.vector(X) == TRUE) | (ncol(X) == 1)) {
        xtrain <- as.vector(Xtrain)
        coef.ll <- lm.ridge.univariate(xtrain, ytrain, lambda = lambda, 
                                       scale = scale)
      }
      else {
        ll <- lm.ridge(ytrain ~ Xtrain, scale = scale, lambda = lambda)
        coef.ll <- coef(ll)
      }
      res <- matrix( length(ytest), length(lambda))
      pred <- t(matrix(coef.ll[, 1], nrow = length(lambda), 
                       ncol = length(ytest))) + Xtest %*% t(coef.ll[, -1])
      res <- pred - matrix(ytest, nrow = length(ytest), ncol = length(lambda))
      cv <- cv + apply(res^2, 2, sum)
    }
    cv <- cv/n
    lambda.opt <- lambda[which.min(cv)]
    if (plot.it == TRUE) {
      plot(lambda, cv, type = "l")
    }
    if ((is.vector(X) == TRUE) | (ncol(X) == 1)) {
      x <- as.vector(X)
      coefficients <- as.vector(lm.ridge.univariate(x, y, scale = scale, 
                                                    lambda = lambda.opt))
    }
    else {
      rr <- lm.ridge(y ~ X, scale = scale, lambda = lambda.opt)
      coefficients <- coef(rr)
    }
    intercept <- coefficients[1]
    coefficients <- coefficients[-1]
    return(list(intercept = intercept, coefficients = coefficients, 
                lambda.opt = lambda.opt))
  }
  ridge12.cv<-function(X,y,lambda1=2,k=10,vvv){
    lambda1<-lambda1
    ss<-seq(-10,-1,length=500)
    ss<-10^ss
    n<-nrow(X)
    nn<-n- floor(n/k)
    lambda<-ss*nn*ncol(X)
    fittt<-function (lambda){
      samplecovmatrix<-t( Xtrain)%*% Xtrain
      samplecovmatrix[-vvv,-vvv]<- samplecovmatrix[-vvv,-vvv]+lambda1[-vvv,-vvv]
      samplecovmatrix[vvv,vvv]<- samplecovmatrix[vvv,vvv]+lambda*(diag(1,length(vvv),length(vvv)))
      ll<-solve(samplecovmatrix)%*%t(Xtrain)%*%ytrain
      coef.ll<-ll
      return(coef.ll)}
    cv<-rep(0,length(lambda))
    n<-nrow(X)
    all.folds <- split(sample(1:n), rep(1:k,length=n))
    for (i in seq(k)) {
      omit <- all.folds[[i]]
      Xtrain=X[-omit,,drop=FALSE]
      ytrain=y[-omit]
      Xtest=X[omit,,drop=FALSE]
      ytest=y[omit]
      samplecovmatrix<-t( Xtrain)%*% Xtrain
      coef.ll<- lapply(lambda, fittt)
      coef.ll <- matrix(unlist(coef.ll), ncol = ncol(X), byrow = TRUE)
      res<-matrix(length(ytest),length(lambda))
      pred<- Xtest%*%t(coef.ll)
      res<-pred-matrix(ytest,nrow=length(ytest),ncol=length(lambda))
      cv<-cv+apply(res^2,2,sum)
    }
    cv<-cv/n
    lambda.opt<-lambda[which.min(cv)]
    return(list(lambda.opt=lambda.opt))
  }
  graphhh<-function(data,kg1, sq1){
    lam1<-seq(sq1[1],sq1[2],sq1[3])
    n <- nrow(data)
    ncols <- ncol(data)
    X<-data[,1:ncols-1]
    loglike<-c()
    aaa<-list()
    for(g in 1:length(lam1)){
      lam<-lam1[g]
      a<-CVglasso(X,lam=lam, K =kg1)
      aa<-a$Omega
      aa[aa!=0]=1
      diag(aa)<-0
      dimnames(aa) <- list(c(seq(1:ncol(X))),c(seq(1:ncol(X))))
      if(isSymmetric(aa)==TRUE) {
        chest.rip<- mpd(aa)
        Separators=chest.rip$separators
        Separators<-Separators[lapply(Separators,length)>0]
        if(length(Separators)>0){
          aaa[[g]]<-aa
          loglike[g]<-a$Loglik }}}
    if(length(which(sapply(aaa, is.null)))>0){
      aaa<-aaa[-which(sapply(aaa, is.null))]
      loglike<- na.omit(loglike)
      mat <- aaa[[which.max(loglike)]]}
    else{
      mat <- aaa[[which.max(loglike)]]}
    return(mat)}
  ridge.group<-function(data,k,PE1,aaa){
    n <- nrow(data)
    ncols <- ncol(data)
    column <- ncols-1
    X <- data[,-ncols]
    y <- as.matrix(data[,ncols])
    chest.rip<- mpd(aaa)
    Cliques=chest.rip$cliques
    C<-list()
    for(l in 1:length(chest.rip$separators)){
      C[[l]] <-  as.numeric(c(do.call("cbind",Cliques[l])) )
      CC <-  as.numeric(c(do.call("cbind",chest.rip$separators[l])) )
      if(length(CC)>0){
        C[[l]]<- C[[l]][-which(C[[l]] %in% CC)]
      }
      else{
        C[[l]]<-C[[l]]}}
    C<-C[lapply(C,length)>0]
    lambda1<-(ridge.cv(X,y,k=k))$lambda.opt
    mse = function(x,y) { mean((x-y)^2)}
    sse1<-c()
    for(ll in 1:length(C)){
      kk<-C[[ll]]
      sse<-c()
      for(i in 1:PE1){
        training <- sample(1:n,round(n*0.75))
        testing <- (1:n)[-training]
        ridge.object<-ridge.cv(X[training,kk],y[training,],k=k)
        B<-ridge.object$coefficients
        y.hat <- as.matrix(X[testing,kk])%*%as.vector(B)
        sse [i]<- sum((y[testing]-y.hat)^2)}
      sse1[ll]<-mean(sse)}
    kkk<-order(sse1)
    CCC<-list()
    for(f in 1:length(kkk) ){
      lg<-kkk[f]
      CCC[[f]]<-C[[lg]]}
    lmd<-c()
    lambda1<-(ridge.cv(X,y,k=k))$lambda.opt
    for(ff in 1:length(kkk)){
      if(ff==1){
        lambda11<-matrix(0,ncol(X),ncol(X))
        diag(lambda11)<-lambda1}
      if(ff>1){
        if(length(CCC[[ff-1]])>1){
          diag(lambda11[CCC[[ff-1]],CCC[[ff-1]]])<-unlist(lmd[ff-1])}
        else{
          lambda11[CCC[[ff-1]],CCC[[ff-1]]]<-unlist(lmd[ff-1])}}
      lmd[ff]<-ridge12.cv(X,y,lambda1=lambda11,k=k,vvv=CCC[[ff]])
    }
    samplecovmatrix<-t(X)%*%X
    for(fff in 1:length(kkk)){
      samplecovmatrix[ CCC[[fff]], CCC[[fff]]]<- samplecovmatrix[ CCC[[fff]], CCC[[fff]]]+lmd[[fff]]*(diag(1,length(CCC[[fff]]),length(CCC[[fff]])))
    }
    BB<-solve(samplecovmatrix)%*%t(X)%*%y
    y.hatBB <- as.matrix(X)%*%as.vector(BB)
    ssesBB <- (sum((y-y.hatBB)^2))/n
    return(list("Coefficients"=BB,"MSE"=ssesBB,"lambda.opt"= unlist(lmd)))
  }
  aa1<-graphhh(data,kg1=kg, sq1=sq)
  result<-ridge.group(data=data,k=k,PE1=PE,aaa= aa1)
  return(result)}
