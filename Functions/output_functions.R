### Output Functions

## General Forms
fun.out <- function(mm,w,tempy,names_pop = NULL,
                    svydesign = NULL,...){
  # mm, w as numeric vectors
  temp <- numeric(15)
  temp[1] <- sum(w*mm)/sum(w) # mean
  temp[2] <- sum(w^2*(mm - temp[1])^2)/(sum(w)^2) # variance
  temp[3] <- temp[1] - qnorm(0.975) * sqrt(temp[2])
  temp[4] <- temp[1] + qnorm(0.975) * sqrt(temp[2])
  temp[5] <- temp[1] - tempy
  temp[6] <- temp[5]^2
  temp[7] <- tempy >= temp[3] & tempy <= temp[4]
  temp[8] <- sum(w)
  temp[9] <- var(w)
  temp[10] <- mean(w<=0)
  if(!is.null(svydesign)){
    temp[11] <- vcov(svymean(~y, svydesign))
    temp[12] <- confint(svymean(~y, svydesign))[1]
    temp[13] <- confint(svymean(~y, svydesign))[2]
    temp[14] <- tempy >= temp[12] & tempy <= temp[13]
    temp[15] <- object.size(svydesign)
  }else {temp[11:15] <- NA}
  if(!is.null(names_pop)){
    names(temp) <- names_pop
  }
  return(temp)
}


## For Suvvey results
svymean.out <- function(svydesign,...){
  # mm, w as numeric vectors
  temp <- numeric(8)
  temp[1] <- coef(svymean(~y, svydesign))
  temp[2] <- vcov(svymean(~y, svydesign))
  temp[3] <- confint(svymean(~y, svydesign))[1]
  temp[4] <- confint(svymean(~y, svydesign))[2]
  temp[5] <- tempy >= temp[3] & tempy <= temp[4]
  temp[6] <- sum(weights(svydesign))
  temp[7] <- var(weights(svydesign))
  temp[8] <- mean(weights(svydesign)<=0)
  return(temp)
}






