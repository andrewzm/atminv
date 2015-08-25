#' @title Line-plot theme
#' 
#' @description Formats a ggplot object for neat plotting.
#'
#' @return Object of class \code{ggplot}
#' @keywords ggplot
#' @export
#' @examples
#' \dontrun{
#' X <- data.frame(x=runif(100),y = runif(100), z = runif(100))
#' LinePlotTheme() + geom_point(data=X,aes(x,y,colour=z))
#'}
LinePlotTheme <- function() {
    g <- ggplot() + theme(panel.background = element_rect(fill='white', colour='black'),text = element_text(size=20),
                          panel.grid.major =  element_line(colour = "light gray", size = 0.05),
                          panel.border  = element_rect(fill=NA, colour='black'),
                          plot.margin=unit(c(5,5,5,0),"mm")) 
    return (g)
}