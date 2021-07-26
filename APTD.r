## --------------------
## --------------------
## --------------------
## This file contains functions for Adaptive Penalized Tensor Decomposition. 
## Most function are only suitable for three dimension tensor. For higher-order tensor, codes should be modified adaptively.
## --------------------
## --------------------
## --------------------




## --------------------
## Generate tensor given d u v w.
## --------------------
generate_tensor<-function(d,u,v,w) {
  if (length(d) == 1){
    u<-t(matrix(as.numeric(u)))
    v<-t(matrix(as.numeric(v)))
    w<-t(matrix(as.numeric(w)))
    make_list <- list('mat' = t(u),'mat2' = t(v),'mat3'= t(w) )
    aux <- ttl(as.tensor(array(1,c(1,1,1))), make_list, ms = c(1,2,3))
    return(d*aux)
  }else{
    temp<-NULL
    for(k0 in 1:length(d)) {
      make_list <- list('mat' = as.matrix(u[k0,]),'mat2' = as.matrix(v[k0,]),'mat3'= as.matrix(w[k0,]) )
      aux <- ttl(as.tensor(array(1,c(1,1,1))), make_list, ms = c(1,2,3))
      if (is.null(temp)){temp<-as.tensor(array(0,dim(aux)))}
      temp <- temp + d[k0]*aux
    }
    return(temp)
  }
}

## --------------------
## Generate arbitrary penalty matrix.
## --------------------
create_D<-function(x, ord = 2, gamma = 1) {
  D <- list()
  if(is.null(nrow(x))) {
    m=1
    n=length(x)
    } 
  else {
    m=nrow(x)
    n=ncol(x)
    }
  for (i in 1:m) {
    D[[i]] <- diag( 1/abs(as.numeric(diff(diag(n),differences=ord) %*% as.matrix(x[i,])))^gamma  ) %*% diff(diag(n),differences=ord)
    }
  return(D)
  }


## --------------------
## Generate sequences on log scale.
## --------------------
lseq <- function (from, to, length, decreasing = FALSE) {
  stopifnot(from > 0)
  out <- 10^(seq(log10(from), log10(to), length.out = length))
  out <- out[order(out, decreasing = decreasing)]; 
  return(out)
  }


## --------------------
## Calculate the mse if forecasting method is GAM, which is used in cross-validation function to select the best tuning parameters.
## --------------------
cal_pmse_gam_adj<-function(k, year, bs = "tp", return_pred = FALSE) {
     u = tempresult[[k]]$u
     v = tempresult[[k]]$v
     w = tempresult[[k]]$w
     w_pred <- matrix(NA,nrow = num_components,ncol = year)

     for (k0 in 1:num_components){
          model1<-gam(y~s(x, bs = bs), data = data.frame(y=w[k0,],x=seq(1,length(w[1,]))))
          w_pred[k0,] = predict(model1, newdata = data.frame(x=seq((length(w[k0,])+1),(length(w[k0,])+year))))
          }

     aux<-0
     for (k0 in 1:num_components){
          lizt <- list('mat1' = t(t(u[k0,])),'mat2' = t(t(v[k0,])),'mat3' =t(t(w_pred[k0,])))
          #ggg<-as.array(g[k0])
          #dim(ggg)<-c(1,1,1)
          aux1<-ttl(as.tensor(array(1,c(1,1,1))), lizt, ms = c(1,2,3))
          aux<-aux + tempresult[[k]]$d[k0]*aux1
     }
     #convert it back 
     aux = aux*sd1[[k]] + as.tensor(array(mean1[[k]],dim(aux))) 
     
     e<-sqrt(fnorm(aux-V1[[k]])^2)
     if(!return_pred)
          return(e)
     if(return_pred)
          return(list(e = e, w_pred = w_pred))
     }


## --------------------
## Calculate the mse if forecasting method is ARIMA, which is used in cross-validation function to select the best tuning parameters.
## --------------------
cal_pmse_arima_adj<-function(k, year, return_pred = FALSE) {
     u = tempresult[[k]]$u
     v = tempresult[[k]]$v
     w = tempresult[[k]]$w
     w_pred<-matrix(NA,nrow = num_components,ncol = year)
     
     for (k0 in 1:num_components){
          w_pred1=rwf(as.matrix(w)[k0,],h=year,drift = TRUE)
          w_pred[k0,]<-as.data.frame(w_pred1)[,1]
          }

          
     aux<-0    
     for (k0 in 1:num_components){
          lizt <- list('mat1' =t(t(u[k0,])),'mat2' = t(t(v[k0,])),'mat3' =t(t(w_pred[k0,])))
          #ggg<-as.array(g[k0])
          #dim(ggg)<-c(1,1,1)
          aux1<-ttl(as.tensor(array(1,c(1,1,1))), lizt, ms = c(1,2,3))
          aux<-aux + tempresult[[k]]$d[k0]*aux1
     }
     
     ####convert back
     aux = aux*sd1[[k]] + as.tensor(array(mean1[[k]],dim(aux))) 
     
     e<-sqrt(fnorm(aux-V1[[k]])^2)
     if(!return_pred)
          return(e)
     if(return_pred)
          return(list(e = e, w_pred = w_pred))
     }
    
## --------------------
## Calculate the mse if forecasting method is Linear extrapolation, which is used in cross-validation function to select the best tuning parameters.
## --------------------
cal_pmse_linearextra_adj<-function(k, year, return_pred = FALSE) {
     u = tempresult[[k]]$u
     v = tempresult[[k]]$v
     w = tempresult[[k]]$w
     w_pred <- matrix(NA,nrow = num_components,ncol = year)

     for (k0 in 1:num_components){
          w_pred[k0,] = approxExtrap(x = seq(1,length(w[1,])), y = w[k0,], xout = seq(length(w[k0,])+1,by=1,length=year))$y
          }

     aux<-0
     for (k0 in 1:num_components){
          lizt <- list('mat1' = t(t(u[k0,])),'mat2' = t(t(v[k0,])),'mat3' =t(t(w_pred[k0,])))
          #ggg<-as.array(g[k0])
          #dim(ggg)<-c(1,1,1)
          aux1<-ttl(as.tensor(array(1,c(1,1,1))), lizt, ms = c(1,2,3))
          aux<-aux + tempresult[[k]]$d[k0]*aux1
     }
     
     ####convert back
     aux = aux*sd1[[k]] + as.tensor(array(mean1[[k]],dim(aux))) 
     
     e<-sqrt(fnorm(aux-V1[[k]])^2)
     if(!return_pred)
          return(e)
     if(return_pred)
          return(list(e = e, w_pred = w_pred))
     }


## --------------------
## Penalized Tensor Decomposition(trend-filtering) for single factor.
## --------------------
PTD_tf<-function(tnsr,c1,c2,c3, Niter=1000, tol = 1e-5,ku,kv,kw,u_init=NULL,v_init=NULL,w_init=NULL) {
  #tuning parameters can be 0
  #trendfiltering
  #penalty matrix D is defined by default  
  if ( is.null(u_init)){
    temp = gen_startvals(1,tnsr = tnsr)
    u_init = temp$u
    v_init = temp$v
    w_init = temp$w
  }
  counter <- 0
  diff <- 10
  cw_u <- u_init
  cw_v <- v_init
  cw_w <- w_init
  u<-t(matrix(as.numeric(u_init)))
  v<-t(matrix(as.numeric(v_init)))
  w<-t(matrix(as.numeric(w_init)))
  
  while(counter < Niter & diff > tol) {
    
    # update w
    lizt <- list('mat2' = u,'mat3' =v)
    w =   ttl(tnsr, lizt, ms = c(1,2)) 
    if (c3!=0){
      w = glmgen::trendfilter( as.vector(w@data), k = kw , family = "gaussian", method = "admm", lambda = c3)
      w = predict(w, lambda = c3)
      w = t(as.matrix(w))
      w = w/norm(w,"F")
    }else
    {
      w = t(as.matrix(w@data))
      w =  w/norm(w,"F")}
    
    # update u
    lizt <- list('mat2' = v,'mat3' =w)
    u =   ttl(tnsr, lizt, ms = c(2,3)) 
    if (c1!=0){
      u = glmgen::trendfilter( as.vector(u@data), k = ku , family = "gaussian", method = "admm", lambda = c1)
      u = predict(u, lambda = c1)
      u = t(as.matrix(u))
      u = u/norm(u,"F")
    }
    else
    {
      u = t(as.matrix(u@data))
      u =  u/norm(u,"F")}
    
    # update v
    lizt <- list('mat2' = u,'mat3' =w)
    v =   ttl(tnsr, lizt, ms = c(1,3)) 
    if(c2!=0){
      v =  glmgen::trendfilter( as.vector(v@data), k = kv , family = "gaussian", method = "admm", lambda = c2)
      v = predict(v, lambda = c2)
      v = t(as.matrix(v))
      v = v/norm(v,"F")}
    else
    {
      v = t(as.matrix(v@data))
      v =  v/norm(v,"F")}
    
    diff <- mean(c(as.numeric(abs(cw_u-u)),as.numeric(abs(cw_v-v)),as.numeric(abs(cw_w-w))))
    cw_u <- u
    cw_v <- v
    cw_w <- w
    counter <- counter + 1
  }
  
  
  lizt <- list('mat1' =u,'mat2' = v,'mat3' =w)
  d = ttl(tnsr, lizt, ms = c(1,2,3))  
  
  
  return(list(u = u, v= v, w=w, d=as.numeric(d@data)))
}


## --------------------
## Penalized Tensor Decomposition(trend-filtering) for multiple factors.
## --------------------
multiple_tf_simple<-function(tnsr,num_components,d_hat=NULL,u_init= NULL,v_init=NULL,w_init=NULL,c1=NULL,c2=NULL,c3=NULL,Niter=1500,tol=1e-05,ku,kv,kw,trace=TRUE) {
  ##specific tuning parameters
  if ( is.null(u_init)){
    temp = gen_startvals(num_components,tnsr = tnsr)
    u_init = temp$u
    v_init = temp$v
    w_init = temp$w
  }
  if ( is.null(c1))
     c1 <- rep(0, num_components)
  if ( is.null(c2))
     c2 <- rep(0, num_components)
  if ( is.null(c3))
     c3 <- rep(0, num_components)

  if(identical(d_hat,rep(0,num_components)) | is.null(d_hat)){ 
    for(j in 1:num_components) {
      d_hat[j] = product(tnsr,c(1,2,3), u_init[j,],v_init[j,],w_init[j,]) 
    }}
  
  Y = as.tensor(array(0,dim(tnsr)))
  
  for( j in 1:num_components){
    lizt <- list('mat' = as.matrix(u_init[j,]),'mat2' = as.matrix(v_init[j,]),'mat3'= as.matrix(w_init[j,]) )
    aux =   ttl(as.tensor(array(1,c(1,1,1))), lizt, ms = c(1,2,3))
    
    Y = Y + d_hat[j]*aux 
  }##  closing for j
  
  ##################333
  counter <- 0
  diff <- 10
  cw_u = u_init
  cw_v = v_init
  cw_w = w_init
  ######
  
  while(diff>tol & counter < Niter)
  {  
    for(j in 1:num_components)
    {
      lizt <- list('mat' = as.matrix(u_init[j,]),'mat2' = as.matrix(v_init[j,]),'mat3'= as.matrix(w_init[j,]) )
      aux =   ttl(as.tensor(array(1,c(1,1,1))), lizt, ms = c(1,2,3))
      
      Y  = Y -  d_hat[j]*aux 
      ##############################################
      
      
      Z =  tnsr-Y   
      
      # if(iter==1)
      # 
      temp = PTD_tf(Z,c1[j],c2[j],c3[j],Niter=Niter,tol=tol,ku=ku,kv=kv,kw=kw,u_init = u_init[j,,drop=FALSE], v_init = v_init[j,,drop=FALSE], w_init = w_init[j,,drop=FALSE])
      
      ##############################################
      u_init[j,] = as.vector(temp$u)
      v_init[j,] = as.vector(temp$v)      
      w_init[j,] = as.vector(temp$w)
      d_hat[j] = product(Z,c(1,2,3), u_init[j,],v_init[j,],w_init[j,]  )
      
      
      lizt <- list('mat' = as.matrix(u_init[j,]),'mat2' = as.matrix(v_init[j,]),'mat3'= as.matrix(w_init[j,]) )
      aux =   ttl(as.tensor(array(1,c(1,1,1))), lizt, ms = c(1,2,3))
      
      
      Y = Y + d_hat[j]*aux 
      
    }
    
    diff <- mean(c(as.numeric(abs(u_init - cw_u)),as.numeric(abs(v_init - cw_v)),as.numeric(abs(w_init - cw_w))))
    if(trace)
      message("Iteration: ", counter, "\t Difference: ", round(diff,5))
    
    cw_u = u_init
    cw_v = v_init
    cw_w = w_init
    counter <- counter + 1
  } 
  
  out <- list(u =u_init, v = v_init, w = w_init, d = d_hat, diff=diff)
  return(out)    
}

## --------------------
## Adaptive Penalized Tensor Decomposition for single factor.
## --------------------
PTD_D<-function(tnsr,c1,c2,c3, Niter=1500, tol = 1e-5,Du,Dv,Dw,u_init,v_init,w_init) {
  #tuning parameters can be 0
  #Arbitrary penalty matrix
  # init must be unpenalized results 
  counter <- 0
  diff <- 10

  u = cw_u = u_init
  v = cw_v = v_init
  w = cw_w = w_init
  
  while(counter < Niter & diff > tol) {
    
    # update w
    lizt <- list('mat2' = u,'mat3' =v)
    w =   ttl(tnsr, lizt, ms = c(1,2)) 
    if (c3!=0){
      w = genlasso( as.vector(w@data),D=Dw)
      w = predict(w, lambda = c3)$fit
      w = t(as.matrix(w))
      w = w/norm(w,"F")
    }else
    {
      lizt <- list('mat2' = u,'mat3' =v)
      w =   ttl(tnsr, lizt, ms = c(1,2)) 
      w = t(as.matrix(w@data))
      w =  w/norm(w,"F")}
    
    # update u
    lizt <- list('mat2' = v,'mat3' =w)
    u =   ttl(tnsr, lizt, ms = c(2,3)) 
    if (c1!=0){
      u = genlasso( as.vector(u@data),D=Du)
      u = predict(u, lambda = c3)$fit
      u = t(as.matrix(u))
      u = u/norm(u,"F")
    }else
    {
      lizt <- list('mat2' = v,'mat3' =w)
      u =   ttl(tnsr, lizt, ms = c(2,3)) 
      u = t(as.matrix(u@data))
      u =  u/norm(u,"F")}
    
    # update v
    lizt <- list('mat2' = u,'mat3' =w)
    v =   ttl(tnsr, lizt, ms = c(1,3))
    if(c2!=0){
      v = genlasso( as.vector(v@data),D=Dv)
      v = predict(v, lambda = c2)$fit
      v = t(as.matrix(v))
      v = v/norm(v,"F")}
    else
      {lizt <- list('mat2' = u,'mat3' =w)
      v =   ttl(tnsr, lizt, ms = c(1,3)) 
      v = t(as.matrix(v@data))
      v =  v/norm(v,"F")}
    
    diff <- mean(c(as.numeric(abs(cw_u-u)),as.numeric(abs(cw_v-v)),as.numeric(abs(cw_w-w))))
    cw_u <- u
    cw_v <- v
    cw_w <- w
    counter <- counter + 1
  }
  
  
  lizt <- list('mat1' =u,'mat2' = v,'mat3' =w)
  d = ttl(tnsr, lizt, ms = c(1,2,3))  
  
  
  return(list(u = u, v= v, w=w,d=as.numeric(d@data)))
}


## --------------------
## Adaptive Penalized Tensor Decomposition for single factor.
## --------------------
multiple_D_adaptive<-function(tnsr,num_components,d_hat=NULL,u_init,v_init,w_init,c1=NULL,c2=NULL,c3=NULL,Niter=15,tol=1e-05,Du,Dv,Dw,trace=TRUE)
  # tuning parameters for each dimension are the same
  # inital values cant be null
  # Dv and Dw should be a list of matrices
  # to save computation time temp , only run 15 times
{
  ##specific tuning parameters
  if (!is.list(Du)){Du=list(Du)}
  if (!is.list(Dv)){Dv=list(Dv)}
  if (!is.list(Dw)){Dw=list(Dw)}
  if (length(Dv) != num_components & length(Dw) != num_components & length(Du) != num_components){return(print("wrong penalty matrix"))}
  if(is.null(d_hat)){ 
    for(j in 1:num_components) {
      d_hat[j] = product(tnsr,c(1,2,3), u_init[j,],v_init[j,],w_init[j,]) 
    }}
  Y = as.tensor(array(0,dim(tnsr)))
  
  for( j in 1:num_components){
    lizt <- list('mat' = as.matrix(u_init[j,]),'mat2' = as.matrix(v_init[j,]),'mat3'= as.matrix(w_init[j,]) )
    aux =   ttl(as.tensor(array(1,c(1,1,1))), lizt, ms = c(1,2,3))
    
    Y = Y + d_hat[j]*aux 
  }##  closing for j
  
  ##################333
  cw_u <- u_init
  cw_v <- v_init
  cw_w <- w_init
  counter <- 0
  diff <- 10
  
  ######
  
  while(diff>tol & counter < Niter)
  {  
    for(j in 1:num_components)
    {
      lizt <- list('mat' = as.matrix(u_init[j,]),'mat2' = as.matrix(v_init[j,]),'mat3'= as.matrix(w_init[j,]) )
      aux =   ttl(as.tensor(array(1,c(1,1,1))), lizt, ms = c(1,2,3))
      
      Y  = Y -  d_hat[j]*aux 
      ##############################################
      
      
      Z =  tnsr-Y   
      
      # if(iter==1)
      # 
      temp = PTD_D(Z,c1[j],c2[j],c3[j],Niter=Niter,tol=tol,Du = Du[[j]],Dv = Dv[[j]],Dw = Dw[[j]],
                   u_init = u_init[j,,drop=FALSE], v_init = v_init[j,,drop=FALSE], w = w_init[j,,drop=FALSE])
      
      ##############################################
      u_init[j,] = as.vector(temp$u)
      v_init[j,] = as.vector(temp$v)      
      w_init[j,] = as.vector(temp$w)
      d_hat[j] = temp$d
      
      
      lizt <- list('mat' = as.matrix(u_init[j,]),'mat2' = as.matrix(v_init[j,]),'mat3'= as.matrix(w_init[j,]) )
      aux =   ttl(as.tensor(array(1,c(1,1,1))), lizt, ms = c(1,2,3))
      
      
      Y = Y + d_hat[j]*aux 
      
    }
    
    diff <- mean(c(as.numeric(abs(u_init - cw_u)),as.numeric(abs(v_init - cw_v)),as.numeric(abs(w_init - cw_w))))
    if(trace)
      message("Iteration: ", counter, "\t Difference: ", round(diff,5))
    
    cw_u <- u_init
    cw_v <- v_init
    cw_w <- w_init
    counter <- counter + 1
  } 
  
  out <- list(u = u_init, v = v_init, w = w_init, d= d_hat)
  return(out)    
}