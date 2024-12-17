### Output Functions

## General Forms
fun.out <- function(mm,w,truey,tempy,names_pop = NULL,
                    svydesign = NULL,...){
  # mm, w as numeric vectors
  temp <- numeric(16)
  temp[1] <- sum(w*mm)/sum(w) # mean
  temp[2] <- sum(w^2*(mm - temp[1])^2)/(sum(w)^2) # variance
  temp[3] <- temp[1] - qnorm(0.975) * sqrt(temp[2])
  temp[4] <- temp[1] + qnorm(0.975) * sqrt(temp[2])
  temp[5] <- temp[1] - tempy
  temp[6] <- temp[5]^2
  temp[7] <- truey >= temp[3] & truey <= temp[4]
  temp[8] <- tempy >= temp[3] & tempy <= temp[4]
  temp[9] <- sum(w)
  temp[10] <- var(w)
  temp[11] <- mean(w<=0)
  if(!is.null(svydesign)){
    temp[12] <- vcov(svymean(~y, svydesign))
    temp[13] <- confint(svymean(~y, svydesign))[1]
    temp[14] <- confint(svymean(~y, svydesign))[2]
    temp[15] <- truey >= temp[13] & truey <= temp[14]
    temp[16] <- object.size(svydesign)
  }else {temp[12:16] <- NA}
  if(!is.null(names_pop)){
    names(temp) <- names_pop
  }
  return(temp)
}


## For Suvvey results
svymean.out <- function(svydesign,truey,...){
  # mm, w as numeric vectors
  temp <- numeric(10)
  temp[1] <- coef(svymean(~y, svydesign))
  temp[2] <- vcov(svymean(~y, svydesign))
  temp[3] <- confint(svymean(~y, svydesign))[1]
  temp[4] <- confint(svymean(~y, svydesign))[2]
  temp[5] <- temp[1]-truey
  temp[6] <- temp[5]^2
  temp[7] <- truey >= temp[3] & truey <= temp[4]
  temp[8] <- sum(weights(svydesign))
  temp[9] <- var(weights(svydesign))
  temp[10] <- mean(weights(svydesign)<=0)
  return(temp)
}


### Subgroup Outcome Function
fun.pre0 <- function(dt,w,y,pre_pop = pre_pop,...){
  temp <- apply(dt,2,function(x){
    index <- x==1
    out <- sum(x[index]*w[index]*y[index])/sum(w[index]*x[index])
    out <- rbind(out,sum((x[index]*w[index]*y[index] - out)^2)/(sum(w[index]*x[index])^2)) # variance
    ci <- c(out[1] - qnorm(0.975) * sqrt(out[2]),out[1] + qnorm(0.975) * sqrt(out[2]))
    out <- c(sum(index),out,ci) 
    return(out)
  })
  temp <- data.frame(temp)
  temp <- rbind(temp, pre_pop >= temp[4,] & pre_pop <= temp[5,])
  rownames(temp) <- c("size","mean","variance","LB","UB","Coverage")
  return((temp))
}


fun.pre <- function(dt,w,y,pre_pop = pre_pop,
                    dt2,pre_all = pre_all,...){
  temp <- apply(dt,2,function(x){
    index <- x==1
    out <- sum(x[index]*w[index]*y[index])/sum(w[index]*x[index])
    out <- rbind(out,sum((x[index]*w[index]*y[index] - out)^2)/(sum(w[index]*x[index])^2)) # variance
    ci <- c(out[1] - qnorm(0.975) * sqrt(out[2]),out[1] + qnorm(0.975) * sqrt(out[2]))
    out <- c(sum(index),out,ci) 
    return(out)
  })
  temp <- data.frame(temp)
  temp <- rbind(temp, pre_pop >= temp[4,] & pre_pop <= temp[5,], fun.ee(dt2,w,y,pre_all))
  rownames(temp) <- c("size","mean","variance","LB","UB","Coverage","EE")
  return((temp))
}

}



