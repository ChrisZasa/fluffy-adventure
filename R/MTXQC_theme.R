# MTXQC ggplot2 theme modification

r_mtxqc_theme <- function(base_size){
	theme_bw(base_size) %+replace%
		theme(#axis.text.y.right = element_text(colour = "#000000", size = rel(0.8)),
			axis.title = element_text(size = rel(1.2), face = "italic"),
			plot.title = element_text(size = rel(1.2), colour = "#000000", face = "bold",
															margin = margin(t = 0, r = 0, b = 6.6, l = 0, unit = "pt")),
		  panel.grid.major = element_line(size = rel(.1), colour = "#000000"),
		  panel.grid.minor = element_line(size = rel(.05), colour = "#000000"),
			legend.title = element_text(size = rel(1.2), face = "bold", colour = "#000000"))
}

# 
# mtxqc_theme <- function(){
#   theme_minimal(base_family = "Roboto Condensed") +
#     theme(plot.title = element_text(size = rel(1.5), face = "bold"),
#       plot.subtitle = element_text(size = rel(1.1)),
#       plot.caption = element_text(color = "#777777", vjust = 0),
#       axis.title = element_text(size = rel(1.2), face = "italic"),
#       panel.grid.major = element_line(size = rel(.1), colour = "#000000"),
#       panel.grid.minor = element_line(size = rel(.05), colour = "#000000"),
#       #legend.position = "none"
#     )
# }

theme_set(r_mtxqc_theme(base_size = 8))