#' Estimate HR from an IPE model
#'
#' @param ipe.input as passed to IPE_BW 
#' @param ipe.output output from IPE_BW
#' @param distr distribiton to use - default weibull
#' @param recensor should recensoring be done - default of 1 = yes  
#'
#' @return
#' @export
#'
#' @examples
MIPE.cox <- function(ipe.input, ipe.output, distr="weibull",recensor=1){
  
  # get param
  if (ipe.output$conv=="YES"){
    phi <- ipe.output$survreg.out$coefficients[2]
    phi.unique <- 1
  } else{
    phi <- mean(ipe.output$trt.param[91:100,])
    phi.unique <- 0
  }
  
  scale <- NA
  weib.hr <- NA
  # if unique then derive HR from weibull model
  if (phi.unique){
    scale <- ipe.output$survreg.out$scale
    weib.hr <- exp(-phi / scale)
  }
  
  # derive latent time
  latent.time <- ipe.input$event.time
  censor<-ipe.input$censor.ind
  
  if (recensor==0){
    tmp.index = which(ipe.input$trt.ind ==0 & ipe.input$switch.ind==1)
    latent.time[tmp.index] = ipe.input$switch.time[tmp.index] + exp(-phi)*(ipe.input$event.time[tmp.index] - ipe.input$switch.time[tmp.index])
    latent.censor.ind <- censor
  }
  
  if (recensor==1){
    tmp.index = which(ipe.input$trt.ind ==0 & ipe.input$switch.ind==1)
    latent.time[tmp.index] = ipe.input$switch.time[tmp.index] + exp(-phi)*(ipe.input$event.time[tmp.index] - ipe.input$switch.time[tmp.index])
    if(phi>=0){
      renloc<-which(censor==0 & ipe.input$trt.ind==0 )
      latent.time[renloc]<-exp(-phi)*(ipe.input$cutofftime[renloc])
      reeloc<-which(censor==1 & ipe.input$trt.ind==0 )
      comp<-(latent.time[reeloc]<=exp(-phi)*ipe.input$cutofftime[reeloc])*1
      latent.time[reeloc]<-(1-comp)*(exp(-phi)*ipe.input$cutofftime[reeloc])+comp*(latent.time[reeloc])
    }else{
      renloc<-which(censor==0 & ipe.input$trt.ind==0 )
      latent.time[renloc]<-(ipe.input$cutofftime[renloc])
      reeloc<-which(censor==1 & ipe.input$trt.ind==0 )
      comp<-(latent.time[reeloc]<=ipe.input$cutofftime[reeloc])*1
      latent.time[reeloc]<-(1-comp)*ipe.input$cutofftime[reeloc]+comp*(latent.time[reeloc])
    }
    latent.censor.ind<-censor
  }
  
  cox.ipe <- coxph(Surv(latent.time,latent.censor.ind)~ipe.input$trt.ind)
  
  rc<- list(ipe.event.time = latent.time, ipe.censor.ind = latent.censor.ind, cox.ipe = cox.ipe, phi = phi, phi.unique = phi.unique, scale = scale, weib.hr = weib.hr)
  return(rc)
}