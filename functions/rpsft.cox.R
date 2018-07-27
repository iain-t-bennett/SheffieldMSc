#' Applies cox model to an RPSFT model 
#' 
#' Adjustes standard errors based on g-estimation test
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
#' @param rpsft.output list returned by RPSFT() or RPSFT.2pass
#' @param rpsft.input A data frame containing specified variable names
#' @param Grho Defines the test used in the grid search (0 = logrank, 1 = wilcoxon)
#' @param use.latent.only Logical. Should comparison by T vs U (FALSE) or V vs U (TRUE)
#'
#' @return
#' @export
#'
#' @examples
RPSFT.cox <- function(rpsft.output, rpsft.input, Grho = 0, use.latent.only = TRUE){
  
  # combine the observed and latent times
  df <- cbind(rpsft.input, transmute(rpsft.output$latent, latent.event.time = event.time, latent.censor.ind = censor.ind, psi = rpsft.output$psi.chosen))
  
  # derive counterfactual times from latent event times and as T vs U 
  if (use.latent.only){
    df <- mutate(df,
                 lat.on  = pmin(t.on * exp(psi), latent.event.time),
                 lat.off = ifelse((latent.event.time - lat.on) > 0, latent.event.time - lat.on, 0),
                 cfact.time       = ifelse(trt.ind == 1, lat.on * exp(-psi) + lat.off,  latent.event.time ), 
                 cfact.censor.ind = latent.censor.ind
    )
  } else {
    df <- mutate(df, 
                 cfact.time       = ifelse(trt.ind == 1, event.time,  latent.event.time ),
                 cfact.censor.ind = ifelse(trt.ind == 1, censor.ind,  latent.censor.ind )
    )                 
  }
  
  # derive corrected hazard ratios
  cox.rpsft <- coxph(Surv(cfact.time, cfact.censor.ind) ~ trt.ind, data = df)
  
  # correct standard error using Z value
  cox.rpsft$var.orig <- cox.rpsft$var
  z0 <- survdiff(Surv(event.time, censor.ind) ~ trt.ind, data = df, rho = Grho)$chisq ^ 0.5
  cox.rpsft$var <- matrix((abs(coef(cox.rpsft)) / z0)^2)
  
  rc <- list(cox.rpsft = cox.rpsft, cfact.time = df$cfact.time, cfact.censor.ind = df$cfact.censor.ind, psi.unique = rpsft.output$psi.unique)
  
  return(rc)
  
}
