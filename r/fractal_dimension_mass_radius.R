# Fractal dimension estimation using the mass–radius method
# Author: Aleksandra Pierzga
# Context: Master's thesis project
# Language: R
# Library: imager

# Load library
library(imager)

# --------------------------------------------------------------------------------
# The function analiza_fraktalna_masa_promien(folder, prog_proc) takes user-defined
# arguments:
#
# folder – path to the directory containing images
#
# prog_proc – binarization threshold in percent
#
# The function computes the fractal dimension of objects based on images stored
# in the specified folder using the mass–radius method. Each image is converted
# to grayscale and binarized. The analysis is based on determining the relationship
# between the logarithm of the radius and the logarithm of the mass
# (pixel density within a circle).
#
# The slope of this log–log relationship is interpreted as the fractal dimension
# of the object.
#
# The function returns a numeric vector of fractal dimension values.
# --------------------------------------------------------------------------------


analiza_fraktalna_masa_promien <- function(folder, prog_proc) {
  
  prog <- prog_proc / 100
  
  # Get file paths from the directory

  pliki <- list.files(folder, pattern = "\\.(png|jpg|jpeg)$", 
                      full.names = TRUE, ignore.case = TRUE)
  
  wyniki <- c()
  # Load files 
  for (sciezka in pliki) {       
    img <- load.image(sciezka)
    
    # Convert to grayscale
    if (spectrum(img) == 3) {
      img <- grayscale(img)
    }
    
    # Binarization
    img <- img > prog
    img_wektor <- as.numeric(img)
    
    # Determine the object center
    wymiary <- dim(img)
    srodek <- colMeans(which(img == 1, arr.ind = TRUE))
    srodek_x <- srodek[1]
    srodek_y <- srodek[2]
    
    # Determine radii and distances from the image center
    promienie <- seq(10, min(wymiary[1:2]) / 2, 10)
    x <- rep(1:wymiary[1], times = wymiary[2])
    y <- rep(1:wymiary[2], each = wymiary[1])
    odleglosc <- sqrt((x - srodek_x)^2 + (y - srodek_y)^2)
    
    # Compute pixel mass
    masa <- numeric(length(promienie))
    i <- 1  
    for (p in promienie) {
      masa[i] <- sum(img_wektor[odleglosc <= p])
      i <- i + 1
    }
    
    # Regression log-log
    model <- lm(log(masa) ~ log(promienie))
    wymiar <- coef(model)[2]
    
    wyniki <- c(wyniki, wymiar)
  }
  
  return(wyniki)
}

# Get directory path and threshold from the user
folder <- readline(prompt = "Enter the path to the folder containing images: ")
prog_proc <- as.numeric(readline(prompt = "Enter the binarization threshold in %: "))

# Function call and display of results
wyniki <- analiza_fraktalna_masa_promien(folder, prog_proc)
print(wyniki)
