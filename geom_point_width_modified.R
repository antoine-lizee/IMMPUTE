## Change in geom_point in order to support width.
#
# from https://groups.google.com/forum/#!topic/ggplot2/c78MAUwLXnI

# # Here it is in action: 
#   
#   # Use the new width parameter
#   ggplot(mtcars, aes(x=mpg, y=wt, colour=factor(cyl))) +
#   geom_point(shape=21, width=3)


# Save the original version in case you want to use it again
GeomPoint_old <- ggplot2:::GeomPoint


# Define the new version
library(proto)
GeomPoint_new <- proto(ggplot2:::Geom, {
  objname <- "point"
  
  draw_groups <- function(., ...) .$draw(...)
  draw <- function(., data, scales, coordinates, na.rm = FALSE, ...) {    
    data <- remove_missing(data, na.rm, 
                           c("x", "y", "size", "shape", "width"), name = "geom_point")
    if (empty(data)) return(zeroGrob())
    
    with(coord_transform(coordinates, data, scales), 
         ggname(.$my_name(), pointsGrob(x, y, size=unit(size, "mm"), pch=shape, 
                                        gp=gpar(col=alpha(colour, alpha), fill = alpha(fill, alpha), fontsize = size * .pt,
                                                lex = width)))
    )
  }
  
  draw_legend <- function(., data, ...) {
    data <- aesdefaults(data, .$default_aes(), list(...))
    
    with(data,
         pointsGrob(0.5, 0.5, size=unit(size, "mm"), pch=shape, 
                    gp=gpar(
                      col=alpha(colour, alpha), 
                      fill=alpha(fill, alpha), 
                      fontsize = size * .pt, lex = width),
         )
    )
  }
  
  icon <- function(.) {
    pos <- seq(0.1, 0.9, length=6)
    pointsGrob(x=pos, y=pos, pch=19, gp=gpar(col="black", cex=0.5), default.units="npc")
  }
  
  default_stat <- function(.) StatIdentity
  required_aes <- c("x", "y")
  default_aes <- function(.) aes(shape=16, colour="black", size=2, fill = NA, alpha = 1, width = 0.6)
  
})


# Make the new function run in the same environment
environment(GeomPoint_new) <- environment(ggplot2:::GeomPoint)

# Replace ggplot2:::GeomPoint with GeomPoint_new
assignInNamespace("GeomPoint", GeomPoint_new, ns="ggplot2")




