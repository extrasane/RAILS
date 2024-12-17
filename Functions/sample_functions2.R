### Sample Functions
sample.function <- function(nsiz,alpha,beta,gamma,...){
  ##### Population Setup #####
  ##### US pop setup
  ## age
  x1 <- rnorm(nsiz, mean = 20, sd = 5)
  # sum(x1>85)/nsiz
  # sum(x1<20)/nsiz
  # plot(density(x1))
  catx1 <- (x1>= quantile(x1,probs = 0.3))*1 + (x1>= quantile(x1,probs = 0.75))*1
  catx1 <- factor(catx1)
  
  ## Sex
  x2 <- rbinom(nsiz,size = 1,prob = 0.65) #1 as female
  catx2 <- factor(x2)
  # sex <- factor(x2,labels = c("Male","Female"))
  
  ## Income in k
  x3 <- rtruncpareto(nsiz, lower = 8,shape = 1,upper = 300)
  qx3 <- quantile(x3,probs = c(0.1,0.5,0.9))
  x3 <- (x3>=qx3[1] & x3<qx3[2])*1 + (x3>=qx3[2] & x3<qx3[3])*2 +
    (x3>=qx3[3])*3
  catx3 <- factor(x3)
  # income <- factor(x3,labels = c("<10k","10k-100k","100k-200k",">200k"))
  
  ### Extra variables
  
  x5 <- sample(c(0,1,2,3),nsiz,replace = TRUE,prob = c(0.29,0.3,0.45,0.005))
  catx5 <- factor(x5)
  x4 <- rbinom(nsiz,size = 1,prob = 0.1*x5)
  catx4 <- factor(x4)
  
  
  ## Occurrence Probability/Prevalence
  # Preset Alpha to determine the variable prevalence
  py <- expit(alpha[1] + alpha[2] * x1 + alpha[3]*x2 + alpha[4]*I(x3==1) +
                alpha[5]*I(x3==2) + alpha[6]*I(x3==3) + alpha[7]*x1*x2 + 
                alpha[8]*x2*I(x3) + alpha[9]*x4 + alpha[10]*I(x5==1) +
                alpha[11]*I(x5==2) + alpha[12]*I(x5==3)
  )
  
  y <- rbinom(nsiz,size = 1, prob = py)
  dt <- data.frame(1:nsiz,catx1,catx2,catx3,catx4,catx5,y,py)
  colnames(dt) <- c("id","catx1","catx2","catx3","catx4","catx5","y","py")
  
  ##### NHANES pop
  ## Preset Beta to determine the sampling probability for AoU pop
  ps <- expit(beta[1] + beta[2] * x1 + beta[3]*x2 + beta[4]*I(x3==1) +
                beta[5]*I(x3==2) + beta[6]*I(x3==3) + beta[7]*x1*x2 + 
                beta[8]*x2*I(x3) + beta[9]*x4 + beta[10]*I(x5==1) +
                beta[11]*I(x5==2) + beta[12]*I(x5==3)
              )
  s <- rbinom(nsiz,size = 1, prob = ps)
  #sum(s)
  dt <- cbind(dt,s = s,ps = ps,fpc = nsiz)
  
  
  ##### AoU 
  ## paou as sampling probability
  ## aou = 1 means included in AoU
  
  paou <- expit(gamma[1] + gamma[2] * x1 + gamma[3]*x2 + gamma[4]*I(x3==1) +
                  gamma[5]*I(x3==2) + gamma[6]*I(x3==3) + gamma[7]*x1*x2 + 
                  gamma[8]*x2*I(x3) + gamma[9]*x4 +  gamma[10]*I(x5==1) +
                  gamma[11]*I(x5==2) + gamma[12]*I(x5==3)

                )
  aou <- rbinom(nsiz,size = 1, prob = paou)
  dt <- cbind(dt,aou = aou,paou)
  # sum(dt$aou==1)
  return(dt)
}           