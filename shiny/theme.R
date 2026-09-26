pal_cat <- c("#00FA9A", "#1F7A8C", "#FFB100", "#D81E5B", "#6A4C93")
pal_seq <- c("#E6FFF5", "#7FFDCD", "#00FA9A", "#00A870", "#00734D", "#04251C")
ink     <- "#00A870"

scale_colour_spring <- function(...) ggplot2::scale_colour_manual(values = pal_cat, ...)
scale_color_spring  <- scale_colour_spring
scale_fill_spring   <- function(...) ggplot2::scale_fill_manual(values = pal_cat, ...)
scale_fill_spring_c <- function(...) ggplot2::scale_fill_gradientn(colours = pal_seq, ...)

# ink, not #00FA9A: spring green is ~1.3:1 on white and unreadable as text
app_theme <- bslib::bs_theme(version = 5, primary = ink, success = "#00FA9A")

ggplot2::theme_set(ggplot2::theme_minimal(base_size = 13))
