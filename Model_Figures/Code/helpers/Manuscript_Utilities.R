# Script defining some simple utility function for teh creation of figures for the immpute manuscript.
#
# Copyright Antoine Lizee 09/2014 antoine.lizee@gmail.com

printGGplot <- function(plot, file, formats = c("png", "tiff", "pdf"), res, ...) {
  for (format in formats) {
    ggsave(filename = paste(file, format, sep = "." ), plot = plot, dpi = res, ...)
  }
}

cmToInches <- function(x) {
  x / 2.54  
}