#' Calculate latent survival time from an RPSFT model
#'
#' rpsft.input must contain the following variables 
#' \itemize{
#'   \item event.time - time to event (should equal t.on + t.off)
#'   \item trt.ind - randomized treatment indicator (0 = placebo, 1 = active)
#'   \item censor.ind - censoring indicator (0 = censored, 1 = event)
#'   \item t.on  - time on active treatment
#'   \item t.off - time off active treatment
#'   \item cutofftime - time of clinical cutoff used for recensoring
#'   }
#' 
#' Will return a data frame with same number of rows as rpsft.input containing following variables
#' \itemize{
#'   \item cstar - recensoring time on latent time scale
#'   \item event.time - latent event time
#'   \item censor.ind - censoring indicator for latent event time (0 = censored, 1 = event)
#'   }
#'   
#' @param rpsft.input A data frame containing specified variable names
#' @param psi A value of psi to use in estimating latent time
#'
#' @return
#' @export
#'
#' @examples
RPSFT.latent <- function(psi, rpsft.input){
  # recensor time
  latent<-data.frame(cstar = pmin(rpsft.input$cutofftime*exp(psi), rpsft.input$cutofftime))
  # counterfactual time
  latent$event.time <- pmin(rpsft.input$t.off + rpsft.input$t.on*exp(psi), latent$cstar )    
  # apply recensoring
  latent$censor.ind <- as.numeric((latent$event.time < latent$cstar) & rpsft.input$censor.ind)
  return(latent)
}
