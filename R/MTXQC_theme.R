# MTXQC ggplot2 theme modification

r_mtxqc_theme <- function(base_size){
	theme_bw(base_size) %+replace%
		theme(axis.text.y.right = element_text(colour = "#000000", size = rel(0.8)),
			axis.title = element_text(colour = "#000000", size = rel(1.2)),
			plot.title = element_text(size = rel(1.2), colour = "#000000",
															margin = margin(t = 0, r = 0, b = 6.6, l = 0, unit = "pt")),
			legend.title = element_text(size = rel(1.2), face = "bold", colour = "#000000"))
}

theme_set(r_mtxqc_theme(base_size = 8))