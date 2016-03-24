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



#' @title ggplot-dual-axis
#' 
#' @description See http://stackoverflow.com/questions/18989001/how-can-i-put-a-transformed-scale-on-the-right-side-of-a-ggplot2
#' @export
ggplot_dual_axis <- function(lhs, rhs, axis.title.y.rhs = "rotate") {
    # 1. Fix the right y-axis label justification
    rhs <- rhs + theme(axis.text.y = element_text(hjust = 0))
    # 2. Rotate the right y-axis label by 270 degrees by default
    if (missing(axis.title.y.rhs) |
        axis.title.y.rhs %in% c("rotate", "rotated")) {
            rhs <- rhs + theme(axis.title.y = element_text(angle = 270))
        }
    # 3a. Use only major grid lines for the left axis
    lhs <- lhs + theme(panel.grid.minor = element_blank())
    # 3b. Use only major grid lines for the right axis
    #     force transparency of the backgrounds to allow grid lines to show
    rhs <- rhs + theme(panel.grid.minor = element_blank(),
                       panel.background = element_rect(fill = "transparent", colour = NA),
                       plot.background = element_rect(fill = "transparent", colour = NA))
    # Process gtable objects
    # 4. Extract gtable
    library("gtable") # loads the grid package
    g1 <- ggplot_gtable(ggplot_build(lhs))
    g2 <- ggplot_gtable(ggplot_build(rhs))
    # 5. Overlap the panel of the rhs plot on that of the lhs plot
    pp <- c(subset(g1$layout, name == "panel", se = t:r))
    g <- gtable_add_grob(g1,
                         g2$grobs[[which(g2$layout$name == "panel")]], pp$t, pp$l, pp$b, pp$l)
    # Tweak axis position and labels
    ia <- which(g2$layout$name == "axis-l")
    ga <- g2$grobs[[ia]]
    ax <- ga$children[["axis"]]  # ga$children[[2]]
    ax$widths <- rev(ax$widths)
    ax$grobs <- rev(ax$grobs)
    ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
    g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1)
    g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)
    g <- gtable_add_grob(g, g2$grobs[[7]], pp$t, length(g$widths), pp$b)
    # Display plot with arrangeGrob wrapper arrangeGrob(g)
    library("gridExtra")
    grid.newpage()
    return(arrangeGrob(g))
}

#' @export
draw_world <- function (g = ggplot(), inc_border = TRUE) 
{
    if (!(is(g, "ggplot"))) 
        stop("g has to be of class ggplot")
    if (!(is.logical(inc_border))) 
        stop("inc_border needs to be TRUE or FALSE")
    long <- lat <- group <- NULL
    data(world, envir = environment(), package = "atminv")
    if (inc_border) {
        border <- data.frame(long = c(-179.99, -179.99, 179.99, 
                                      179.99), lat = c(-89.99, 89.99, 89.99, -89.99), group = 1e+05, 
                             region = "border")
        world <- plyr::rbind.fill(world, border)
    }
    g + geom_path(data = world, aes(x = long, y = lat, group = group), 
                  colour = "black", size = 0.1)
}
