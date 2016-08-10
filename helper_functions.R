waterfall_pst <-function(X, k=NULL, seed=1, x.reverse=F){
  y <- X
  # r <- prcomp(t(x))
  # y <- r$x*matrix(r$sdev^2/sum(r$sdev^2),nrow=nrow(r$x),ncol=ncol(r$x),byrow=T)
  #y <-y[order(y[,1]),]
  #u <- r$rotation
  
  # kmeans
  set.seed(seed)
  r <- kmeans(y,k)
  z <- r$centers
  z <- z[order(z[,1]),]
  rownames(z) <-paste0("t",1:nrow(z))
  m <- ape::mst(dist(z))
  
  t.names <-names(which(colSums(m!=0)==1))[1] # There are two ends, then use the left most one.
  for (i in 1:nrow(m)){
    t.names <-append(t.names,names(which(m[t.names[i],]==1))[which(names(which(m[t.names[i],]==1)) %notin% t.names)])
  }
  
  y2d <-y[,1:2]
  #y2d <-y2d[order(y2d[,1]),]
  z2d <-z[,1:2]
  z2d <-z2d[t.names,]
  
  time_start.i <-0
  updatethis.dist <-rep(Inf,nrow(y2d))
  updatethis.time <-rep(0,nrow(y2d))
  update.updown <-rep(0,nrow(y2d))
  pseudotime.flow <-c(0)
  
  for (i in 1:(nrow(z2d)-1)){
    
    # distance between this z2d.i and all y2d
    dot.dist.i <-apply(y2d,1,function(X){dist(rbind(X,z2d[i,]))})
    
    # distance between this z2d.i-z2d.i+1 segment and "insider" y2d
    inside_this_segment <-which(apply(y2d,1,function(X){inside_check.foo(z2d[i,],z2d[i+1,],X)}))
    seg.dist.i <-rep(Inf,nrow(y2d))
    seg.dist.i[inside_this_segment] <-apply(y2d,1,function(X){distance.foo(z2d[i,],z2d[i+1,],X)})[inside_this_segment]
    
    # intersect coordinate between this z2d.i-z2d.i+1 segment and all y2d
    intersect.i <-t(apply(y2d,1,function(X){intersect.foo(z2d[i,],z2d[i+1,],X)}))
    
    # this z2d.i-z2d.i+1 segment's unit vector
    seg_unit_vector <-unit_vector.foo(z2d[i,],z2d[i+1,])
    
    # UPDATE
    # 2. idx for the shortest distance at this round (either dot or seg)
    update.idx <-apply(cbind(dot.dist.i,seg.dist.i,updatethis.dist),1,which.min)
    # 3. update the pseudotime for y2ds with the short distance from the z2d.i
    updatethis.time[which(update.idx==1)] <-time_start.i
    # 4. update the pseudotime for y2ds with the short distance from the z2d.i-z2d.i+1 segment
    relative_cordinates <-t(apply(intersect.i[which(update.idx==2),],1,function(X){seg_unit_vector%*%(X-z2d[i,])}))
    updatethis.time[which(update.idx==2)] <-time_start.i + relative_cordinates
    # 1. update the shortest distance
    updatethis.dist <-apply(cbind(dot.dist.i,seg.dist.i,updatethis.dist),1,min)
    
    update.updown[which(update.idx==1)] <-c(apply(y2d,1,function(X){crossvec_direction(z2d[i,],z2d[i+1,],X)})*dot.dist.i)[which(update.idx==1)]
    update.updown[which(update.idx==2)] <-c(apply(y2d,1,function(X){crossvec_direction(z2d[i,],z2d[i+1,],X)})*seg.dist.i)[which(update.idx==2)]
    
    # update time for the next round
    time_start.i <-time_start.i + dist(rbind(z2d[i,],z2d[i+1,]))
    pseudotime.flow <-append(pseudotime.flow,time_start.i)
  }
  
  # For the y2ds that are closest to the starting z2d
  i=1
  dot.dist.i <-apply(y2d,1,function(X){dist(rbind(X,z2d[i,]))})
  if (length(start.idx <-which(dot.dist.i <= updatethis.dist))>0){
    intersect.i <-t(apply(y2d,1,function(X){intersect.foo(z2d[i,],z2d[i+1,],X)}))
    seg_unit_vector <-unit_vector.foo(z2d[i,],z2d[i+1,])
    relative_cordinates <-0 + t(apply(intersect.i,1,function(X){seg_unit_vector %*% (X-z2d[i,])}))[start.idx]
    updatethis.time[start.idx] <-relative_cordinates
    seg.dist.i <-apply(y2d,1,function(X){distance.foo(z2d[i,],z2d[i+1,],X)})
    update.updown[start.idx] <-c(apply(y2d,1,function(X){crossvec_direction(z2d[i,],z2d[i+1,],X)})*seg.dist.i)[start.idx]
  }
  # For the y2ds that are closest to the arriving z2d
  i=nrow(z2d)
  dot.dist.i <-apply(y2d,1,function(X){dist(rbind(X,z2d[i,]))})
  if (length(arrive.idx <-which(dot.dist.i <= updatethis.dist))>0){
    intersect.i <-t(apply(y2d,1,function(X){intersect.foo(z2d[i-1,],z2d[i,],X)}))
    seg_unit_vector <-unit_vector.foo(z2d[i-1,],z2d[i,])
    relative_cordinates <-time_start.i + as.numeric(t(apply(intersect.i,1,function(X){seg_unit_vector %*% (X-z2d[i,])})))[arrive.idx]
    updatethis.time[arrive.idx] <-relative_cordinates
    seg.dist.i <-apply(y2d,1,function(X){distance.foo(z2d[i-1,],z2d[i,],X)})
    update.updown[arrive.idx] <-c(apply(y2d,1,function(X){crossvec_direction(z2d[i-1,],z2d[i,],X)})*seg.dist.i)[arrive.idx]
  }
  
  pseudotime <-updatethis.time
  pseudotime.y <-update.updown
  pseudotime.flow <-pseudotime.flow
  
  if (x.reverse==T){
    pseudotime <- -pseudotime
    pseudotime.flow <- -pseudotime.flow
  }
  
  pseudotime_range <-max(pseudotime)-min(pseudotime)
  
  pseudotime.flow <-pseudotime.flow-min(pseudotime)
  pseudotime.flow <-pseudotime.flow/pseudotime_range
  
  pseudotime <-pseudotime-min(pseudotime)
  pseudotime <-pseudotime/pseudotime_range

  return(pseudotime)
  
}

slingshot_pst <- function(X, clus.labels, ...){
  l <- get_lineages(X, clus.labels, ...)
  c <- get_curves(X, clus.labels, l, ...)
  out <- sapply(c,function(i){i$pseudotime})
  return(out)
}

monocle_pst <- function(X, reverse=F){
  require(monocle)
  
}

# s_{\pi_1\pi_2} from TSCAN paper
Spp <- function(pst1, pst2){
  pst1 <- as.numeric(pst1)
  pst2 <- as.numeric(pst2)
  # nicely formatted pst1, pst2
  spp <- function(p1,p2){
    A <- length(p1)
    cons <- 2/(A*(A-1))
    sigma <- sapply(seq_len(A),function(i){
      sapply(seq_len(A-i),function(k){
        j <- i+k
        if(any(is.na(c(p1[c(i,j)],p2[c(i,j)])))){
          return(0)
        }
        if((p1[j]-p1[i])*(p2[j]-p2[i]) > 0){
          return(1)
        }
        return(0)
      })
    })
    return(cons*sum(unlist(sigma)))
  }
  if(is.null(names(pst1)) | is.null(names(pst2))){
    if(length(pst1)==length(pst2)){
      s <- spp(pst1,pst2)
    }else{
      stop('lengths must match or names must be provided')
    }
  }else{
    ns <- union(names(pst1),names(pst2))
    p1 <- pst1[match(ns,names(pst1))]
    p2 <- pst2[match(ns,names(pst2))]
    s <- spp(p1,p2)
  }
  return(s)
}