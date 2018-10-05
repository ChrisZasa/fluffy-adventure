#### MTXQC color codes

### some inspiration: https://www.canva.com/learn/website-color-schemes/


library(ggplot2)

theme_set(theme_bw())


# drsimonj corporate colors
mtxqc_colors <- c(
  `orangecream`        = "#F69256",
  `orange` = "#FDBF54",
  `apricotlight` = "#FEE6B2",
  `jade`       = "#6cc8be",
  `icedmint`     = "#bfe2d2",
  `emarald`     = "#007E6F",
  `seafoam` = "#6ACCE8",
  `pear`  = "#90CA79",
  `turquoise` = "#1CB89C",
  `fire` = "#DC4D25",
  `darkred` = "#A61D21",
  `teal` = "#004A70",
  `blue` = "#038589", 
  `darkblue` = "#004A70", 
  `ligthgrey` = "#CACACB",
  `darkgrey` = "#858586",
  `navyblue` = "#18121E",
  `gunmetal` = "#233237",
  `rustyred` = "#984B43",
  `warm_yellow` = "#EAC67A",
  `ice` = "#99D3DF",
  `freshwater` = "#88BBD6",
  `plaster` = "#CDCDCD",
  `linen` = "#E9E9E9",
  `fresh` = "#F7F5E6",
  `vermillion` = "#333A56",
  `sunshine` = "#52658F",
  `clean` = "#E8E8E8", 
  `blueberry` = "#6B7A8F", 
  `apricot` = "#F7882F", 
  `citrus` = "#F7C331", 
  `applecore` = "#DCC7AA",
  `sunflowers` = "#FECE00", 
  `starrynight` = "#0375B4", 
  `irises` = "#007849", 
  `evening` = "#262228"
  
)

#' Function to extract drsimonj colors as hex codes
#'
#' @param ... Character names of drsimonj_colors 
#'
mtxqc_cols <- function(...) {
  cols <- c(...)
  
  if (is.null(cols))
    return(mtxqc_colors)
  
  mtxqc_colors[cols]
}


##### Create palettes

mtxqc_palettes <- list(
  `main`  = mtxqc_cols("orangecream", "orange", "apricotlight", "turquoise", "lightgrey", "darkgrey"),
  
  `cool`  = mtxqc_cols("orangecream", "apricot", "jade", "icedmint"),
  
  `hot`   = mtxqc_cols("blue", "seafoam", "apricotlight", "fire", "darkred", "darkgrey") ,
  
  `akaru` = mtxqc_cols("warm_yellow", "rustyred", "gunmetal", "navyblue"), 
  
  `clean` = mtxqc_cols("ice", "freshwater", "plaster", "linen"),
  
  `sleek` = mtxqc_cols("clean", "sunshine", "vermillion"),
  
  `greens` = mtxqc_cols("emarald", "seafoam", "pear"),
  `cheerful` = mtxqc_cols("blueberry", "apricot", "citrus", "applecore"), 
  `art` = mtxqc_cols("sunflowers", "starrynight", "irises", "evening")
  
)


#' Return function to interpolate a drsimonj color palette
#'
#' @param palette Character name of palette in drsimonj_palettes
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments to pass to colorRampPalette()
#'
mtxqc_pal <- function(palette = "main", reverse = FALSE, ...) {
  pal <- mtxqc_palettes[[palette]]
  
  if (reverse) pal <- rev(pal)
  
  colorRampPalette(pal, ...)
}


##### Create scales

#' Color scale constructor for drsimonj colors
#'
#' @param palette Character name of palette in drsimonj_palettes
#' @param discrete Boolean indicating whether color aesthetic is discrete or not
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_color_gradientn(), used respectively when discrete is TRUE or FALSE
#'
scale_color_mtxqc <- function(palette = "main", discrete = TRUE, reverse = FALSE, ...) {
  pal <- mtxqc_pal(palette = palette, reverse = reverse)
  
  if (discrete) {
    discrete_scale("colour", paste0("mtxqc_", palette), palette = pal, ...)
  } else {
    scale_color_gradientn(colours = pal(256), ...)
  }
}


#' Fill scale constructor for drsimonj colors
#'
#' @param palette Character name of palette in drsimonj_palettes
#' @param discrete Boolean indicating whether color aesthetic is discrete or not
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_fill_gradientn(), used respectively when discrete is TRUE or FALSE
#'
scale_fill_mtxqc <- function(palette = "main", discrete = TRUE, reverse = FALSE, ...) {
  pal <- mtxqc_pal(palette = palette, reverse = reverse)
  
  if (discrete) {
    discrete_scale("fill", paste0("mtxqc_", palette), palette = pal, ...)
  } else {
    scale_fill_gradientn(colours = pal(256), ...)
  }
}


ggplot(iris, aes(Sepal.Width, Sepal.Length, color = Sepal.Length)) +
  geom_point(size = 4) +
  scale_color_mtxqc(discrete = FALSE, palette = "main")

ggplot(iris, aes(Sepal.Width, Sepal.Length, color = Sepal.Length)) +
  geom_point(size = 4) +
  scale_color_mtxqc(discrete = FALSE, palette = "clean")

ggplot(iris, aes(Sepal.Width, Sepal.Length, color = Sepal.Length)) +
  geom_point(size = 4) +
  scale_color_mtxqc(discrete = FALSE, palette = "cool")

ggplot(iris, aes(Sepal.Width, Sepal.Length, color = Sepal.Length)) +
  geom_point(size = 4) +
  scale_color_mtxqc(discrete = FALSE, palette = "greens")

ggplot(iris, aes(Sepal.Width, Sepal.Length, color = Sepal.Length)) +
  geom_point(size = 4) +
  scale_color_mtxqc(discrete = FALSE, palette = "akaru")

ggplot(iris, aes(Sepal.Width, Sepal.Length, color = Sepal.Length)) +
  geom_point(size = 4) +
  scale_color_mtxqc(discrete = FALSE, palette = "sleek", reverse = TRUE)



ggplot(mpg, aes(manufacturer, fill = manufacturer)) +
  geom_bar(size = .2, color = "black") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_mtxqc(palette = "sleek", guide = "none")

ggplot(mtcars, aes(hp, mpg)) +
  geom_point(color = mtxqc_cols("emarald"),
             size = 4, alpha = .8)
