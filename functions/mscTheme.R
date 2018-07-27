
#' Theme for MSc plots
#'
#' @param base_size default - 24
#' @param ... passed to theme_classic()
#'
#' @return ggplot2 theme
#' @export ggplot2
#' 
mscTheme <- function(base_size = 24, ...){
  
  
  rc <-  theme_classic(base_size = base_size, ...) +
    theme(
      axis.line.y = element_line(size = 1, colour = "black"),
      axis.line.x = element_line(size = 1, colour = "black"),
      axis.title = element_text( colour = "black",face = "bold"),
      axis.text.x = element_text( colour = "black",face = "bold"),
      axis.text.y = element_text( colour = "black",face = "bold"),
      legend.text = element_text( colour = "black",face = "bold"),
      axis.ticks.length = unit(0,"cm"),
      plot.margin = unit(c(0,0,0,0), "cm"),
      strip.text.x = element_text(colour = "black",face = "bold"),
      strip.background = element_blank(),
      legend.title = element_blank() 
    ) 

  return(rc)
}