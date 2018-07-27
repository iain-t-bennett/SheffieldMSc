#' Compare latent survival time from an RPSFT model
#'
#' Compares latent survival time from an RPSFT model for a given value of psi
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
#' @param rpsft.input A data frame containing specified variable names
#' @param psi A value of psi to use in estimating latent time
#' @param Grho Defines the test used (0 = logrank, 1 = wilcoxon)
#'
#' @return Z value for the given value of psi
#' @export
#'
#' @examples
RPSFT.trypsi <- function(psi, Grho, rpsft.input){
  
  # define latent survival time given psi
  latent <- RPSFT.latent(psi, rpsft.input)
  
  # derive Z value for chosen test
  sv <- Surv(latent$event.time, latent$censor.ind)
  sd <- survdiff(sv ~ rpsft.input$trt.ind, rho = Grho)
  rc <- sign(sd$obs[2]-sd$exp[2])*sd$chisq^0.5
  return(rc)
}

