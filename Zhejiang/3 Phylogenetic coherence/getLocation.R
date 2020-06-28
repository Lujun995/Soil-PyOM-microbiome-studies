getLocation<-function (p, data, offset = 0, width = 0.1, low = "green", high = "red", 
                                 border_color = NULL, colnames = TRUE, colnames_position = "bottom", 
                                 colnames_angle = 0, colnames_level = NULL, colnames_offset_x = 0, 
                                 colnames_offset_y = 0, font.size = 4, hjust = 0.5) 
{
  #the function is modified from ggtree::gheatmap
  require(magrittr)
  require(tidyr)
  colnames_position %<>% match.arg(c("bottom", "top"))
  variable <- value <- lab <- y <- NULL
  width <- width * (p$data$x %>% range %>% diff)/ncol(data)
  isTip <- x <- y <- variable <- value <- from <- to <- NULL
  df <- p$data
  df <- df[df$isTip, ]
  start <- max(df$x) + offset
  dd <- as.data.frame(data)
  lab <- df$label[order(df$y)]
  dd <- dd[lab, , drop = FALSE]
  dd$y <- sort(df$y)
  dd$lab <- lab
  dd <- gather(dd, variable, value, -c(lab, y))
  i <- which(dd$value == "")
  if (length(i) > 0) {
    dd$value[i] <- NA
  }
  if (is.null(colnames_level)) {
    dd$variable <- factor(dd$variable, levels = colnames(data))
  }
  else {
    dd$variable <- factor(dd$variable, levels = colnames_level)
  }
  V2 <- start + as.numeric(dd$variable) * width
  mapping <- data.frame(from = dd$variable, to = V2)
  mapping <- unique(mapping)
  dd$x <- V2
  dd$width <- width
  dd
}