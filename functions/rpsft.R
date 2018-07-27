#' Does a grid searche to estimate an RPSFT model
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
#' By doing two grid searches is hopeffuly faster than trying all parameter space with fine grid
#' @param rpsft.input A data frame containing specified variable names
#' @param psimat A matrix of values to try in grid search.
#' @param Grho Defines the test used in the grid search (0 = logrank, 1 = wilcoxon)
#'
#' @return
#' @export
#'
#' @examples
RPSFT<-function(rpsft.input,
                    psimat      = matrix(ncol = 1, 
                                         data = seq(from=-2, to=2, by=0.01)),
                    Grho        = 0 
                    
){  
  
  # derive z values for a range of psi
  z <- apply(psimat,c(1),RPSFT.trypsi, Grho = Grho, rpsft.input = rpsft.input)
  
  # chosen value for psi is that which has rank test Z value of 0
  # select based on when sign changes
  x <- sign(z)
  x1 <- x[-1]
  x2 <- x[-length(x)]
  
  # selected psi index (either is 0 or is between a sign change)
  indx         <- which(z==0)
  indx_chg_low <- which(x1!= x2 & x1!=0 & x2!=0)
  indx_chg_hgh <- indx_chg_low+1
  psi.found    <- c(psimat[indx,1], 0.5*(psimat[indx_chg_low,1]+psimat[indx_chg_hgh,1]))
  
  # is it a unique solution?
  psi.unique <- as.numeric(length(psi.found)==1)
  
  # apply proposal of White to take weighted average if is multiple
  psi.chosen <- sum(c(1,rep(c(-1,1),0.5*(length(psi.found)-1)))*psi.found)
  
  # get latent survival time for chosen value
  latent <- RPSFT.latent(psi.chosen, rpsft.input)
  
  # return the tested values, associated Z values, chosen psi and latent survival
  rc <- list(psi.tried  = as.numeric(psimat),
             z          = z,
             psi.found  = psi.found,
             psi.unique = psi.unique,
             psi.chosen = psi.chosen,
             latent     = latent
             )
  return(rc)
}

