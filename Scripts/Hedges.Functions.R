

###### Functions to calculate Hedges effect size for interactions between stressors

###### Modified from Hale et al. (2017) Ecology and Evolution ... 
###### ... "Describing and understanding behavioral responses to multiple ...
###### ... stressors and multiple stimuli"

###### Created for the paper, Orr et al. (2020) XXXX ... 
###### ... "Rapid Evolution Generates Synergism Between ...
###### ... Multiple Stressors and Complicates Ecological Restoration"

###### James Orr
###### 2020/03/08



######################### When there are two stressors ############################### 

hedges.2 <- function(Con, A, B, AB){
  ### Sample Sizes ###
  ncon <- length(Con)
  na <- length(A)
  nb <- length(B)
  nab <- length(AB)
  ### Means ###
  mcon <- mean(Con)
  ma <- mean(A)
  mb <- mean(B)
  mab <- mean(AB)
  ### Standard Deviations ###
  scon <- sd(Con)
  sa <- sd(A)
  sb <- sd(B)
  sab <- sd(AB)
  ###Calculate pooled standard deviations###
  pooledAC <- sqrt(((na-1)*sa^2+(ncon-1)*scon^2)/(na+ncon-2))
  pooledBC <- sqrt(((nb-1)*sb^2+(ncon-1)*scon^2)/(nb+ncon-2))
  pooledAll <- sqrt(((na-1)*(sa^2)+((nb-1)*(sb^2))+
                       ((ncon-1)*(scon^2))+((nab-1)*(sab^2)))/(na+nb+ncon+nab-4))
  ###Calculate J(m) small-sample bias correction###
  Jma <- 1-(3)/((4*(na+ncon-2))-1)
  Jmb <- 1-(3)/((4*(nb+ncon-2))-1)
  Jmi <- 1-(3)/((4*(ncon+na+nb+nab-4))-1)
  ###Calculate individual effects and their sampling variances###
  da <- ((ma-mcon)/pooledAll)*Jma
  db <- ((mb-mcon)/pooledAll)*Jmb
  sv_a <- ((1/na)+(1/ncon)+((da^2)/(2*(na+ncon))))
  sv_b <- ((1/nb)+(1/ncon)+((db^2)/(2*(nb+ncon))))
  ###Calculate 95% confidence intervals (Z value = 1.96) for individual effects 
  CI_a_low <- da-(1.96*sqrt(sv_a))
  CI_a_high <- da+(1.96*sqrt(sv_a))
  CI_b_low <- db-(1.96*sqrt(sv_b))
  CI_b_high <- db+(1.96*sqrt(sv_b))
  ###Calculate interaction effect size and sampling variance
  dI <- (((mab-mcon)-(ma-mcon)-(mb-mcon))*Jmi/(pooledAll)) 
  vI <- ((1/na)+(1/nb)+(1/ncon)+(1/nab)+((dI^2)/(2*(na+nb+ncon+nab))))
  ###Calculate interaction confidence interval
  ICI_low <- dI-(1.96*(sqrt(vI)))
  ICI_high <- dI+(1.96*(sqrt(vI)))
  return(list(da, CI_a_low, CI_a_high, db, CI_b_low, CI_b_high, dI, ICI_low, ICI_high))
}


######################### When there are three stressors ############################### 

hedges.3 <- function(Con, A, B, C, ABC){
  ### Sample Sizes ###
  ncon <- length(Con)
  na <- length(A)
  nb <- length(B)
  nc <- length(C)
  nabc <- length(ABC)
  ### Means ###
  mcon <- mean(Con)
  ma <- mean(A)
  mb <- mean(B)
  mc <- mean(C)
  mabc <- mean(ABC)
  ### Standard Deviations ###
  scon <- sd(Con)
  sa <- sd(A)
  sb <- sd(B)
  sc <- sd(C)
  sabc <- sd(ABC)
  ###Calculate pooled standard deviations###
  pooledAC <- sqrt(((na-1)*sa^2+(ncon-1)*scon^2)/(na+ncon-2))
  pooledBC <- sqrt(((nb-1)*sb^2+(ncon-1)*scon^2)/(nb+ncon-2))
  pooledCC <- sqrt(((nc-1)*sc^2+(ncon-1)*scon^2)/(nc+ncon-2))
  pooledAll <- sqrt(((na-1)*(sa^2)+((nb-1)*(sb^2))+((nc-1)*(sc^2))+
                       ((ncon-1)*(scon^2))+((nabc-1)*(sabc^2)))/(na+nb+nc+ncon+nabc-5))
  ###Calculate J(m) small-sample bias correction###
  Jma <- 1-(3)/((4*(na+ncon-2))-1)
  Jmb <- 1-(3)/((4*(nb+ncon-2))-1)
  Jmc <- 1-(3)/((4*(nc+ncon-2))-1)
  Jmi <- 1-(4)/((5*(ncon+na+nb+nc+nabc-5))-1)
  ###Calculate individual effects and their sampling variances###
  da <- ((ma-mcon)/pooledAll)*Jma
  db <- ((mb-mcon)/pooledAll)*Jmb
  dc <- ((mc-mcon)/pooledAll)*Jmc
  sv_a <- ((1/na)+(1/ncon)+((da^2)/(2*(na+ncon))))
  sv_b <- ((1/nb)+(1/ncon)+((db^2)/(2*(nb+ncon))))
  sv_c <- ((1/nc)+(1/ncon)+((dc^2)/(2*(nc+ncon))))
  ###Calculate 95% confidence intervals (Z value = 1.96) for individual effects 
  CI_a_low <- da-(1.96*sqrt(sv_a))
  CI_a_high <- da+(1.96*sqrt(sv_a))
  CI_b_low <- db-(1.96*sqrt(sv_b))
  CI_b_high <- db+(1.96*sqrt(sv_b))
  CI_c_low <- dc-(1.96*sqrt(sv_c))
  CI_c_high <- dc+(1.96*sqrt(sv_c))
  ###Calculate interaction effect size and sampling variance      
  dI <- (((mabc-mcon)-(ma-mcon)-(mb-mcon)-(mc-mcon))*Jmi/(pooledAll))       
  vI <- ((1/na)+(1/nb)+(1/nc)+(1/ncon)+(1/nabc)+((dI^2)/(2*(na+nb+nc+ncon+nabc))))
  ###Calculate interaction confidence interval
  ICI_low <- dI-(1.96*(sqrt(vI)))
  ICI_high <- dI+(1.96*(sqrt(vI)))
  return(list(dI, ICI_low, ICI_high))
}



