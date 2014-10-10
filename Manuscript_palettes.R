## Palettes used for the different figures of the manuscript
#
# Copyright Antoine Lizee @ UCSF 08/2014 antoine.lizee@ucsf.edu

p1 <- RColorBrewer::brewer.pal(n = 4, name = "Dark2")
p2 <- RColorBrewer::brewer.pal(n = 4, name = "Paired")
colors_perso <-  c(p2[1], p1[2], p2[3], p1[4])
palette_perso <- scale_color_manual(values = colors_perso)
palette_perso_fill <- scale_fill_manual(values = colors_perso)
## Darken it for alphization down the road.
palette_perso_fill2 <- scale_fill_manual(values = adjustcolor(colors_perso, r=0.95, g=0.95, b=0.95 ))
palette_perso2 <- scale_color_manual(values = adjustcolor(colors_perso, r=0.95, g=0.95, b=0.95 ))
