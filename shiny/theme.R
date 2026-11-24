pal_cat <- c("#00FA9A", "#1F7A8C", "#FFB100", "#D81E5B", "#6A4C93")
pal_seq <- c("#E6FFF5", "#7FFDCD", "#00FA9A", "#00A870", "#00734D", "#04251C")
# Subcluster colours: pal_cat, then colours picked far apart by CIEDE2000
pal_wide <- c(pal_cat, "#85660D", "#F99379", "#FE00FA", "#008856", "#90AD1C", "#3283FE",
              "#882D17", "#D55E00", "#E68FAC", "#3B00FB", "#E41A1C", "#325A9B", "#B10DA1",
              "#1CBE4F", "#AA0DFE", "#F28E2B", "#7570B3", "#E15759", "#B00068", "#A65628",
              "#C075A6")
ink     <- "#00A870"

scale_colour_spring <- function(...) ggplot2::scale_colour_manual(values = pal_cat, ...)
scale_color_spring  <- scale_colour_spring
scale_fill_spring   <- function(...) ggplot2::scale_fill_manual(values = pal_cat, ...)
scale_fill_spring_c <- function(...) ggplot2::scale_fill_gradientn(colours = pal_seq, ...)

# The app's Bootstrap theme, ink as primary and a root font size of 85%
app_theme <- bslib::bs_theme(version = 5, primary = ink, success = "#00FA9A",
                             "font-size-root" = "85%")

# The plots' base text size, scaled by the app's figure size control
PLOT_PT <- 13
ggplot2::theme_set(ggplot2::theme_minimal(base_size = PLOT_PT))
