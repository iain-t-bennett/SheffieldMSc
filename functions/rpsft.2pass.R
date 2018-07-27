#' Does 2 grid searches to estimate an RPSFT model
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
#' @param pass1.psimat A matrix of values to try in first grid search. Second will be narrower
#' @param Grho Defines the test used in the grid search (0 = logrank, 1 = wilcoxon)
#'
#' @return
#' @export
#'
#' @examples
RPSFT.2pass <- function(rpsft.input, pass1.psimat = matrix(ncol = 1, data = c(seq(from=-3, to=3, by=0.1))), Grho = 0) {
  
  # apply RPSFT with wide grid
  pass1 <- RPSFT(rpsft.input, Grho=Grho, psimat = pass1.psimat)
  
  # get interesting places to look (within 95 percent CI for z)
  if (any(abs(pass1$z) <= 1.96)){
    pass1.limits <- pass1$psi.tried[abs(pass1$z) <= 1.96]
    pass2.psimat <- matrix(ncol = 1, data = c(seq(from = min(pass1.limits), to = max(pass1.limits), by = 0.01)))
  } else{
    pass2.psimat <- matrix(ncol = 1, data = c(seq(from = -5, to = 5, by = 0.05)))
  }
  
  # apply RPSFT again
  pass2 <- RPSFT(rpsft.input, Grho=Grho, psimat = pass2.psimat)
  
  # combine the grids searched for return
  pass.df <- rbind(data.frame(pass = 1, psi=pass1$psi.tried, z = pass1$z),
                   data.frame(pass = 2, psi=pass2$psi.tried, z = pass2$z))
  
  pass.df <- summarise(group_by(pass.df, psi , z), pass = min(pass))
  
  # return the second pass but add in the values tried in the first pass
  rc <- pass2
  rc$psi.tried <- pass.df$psi
  rc$z         <- pass.df$z
  rc$pass      <- pass.df$pass
  
  return(rc)
}
