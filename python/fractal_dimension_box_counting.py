"""
Fractal dimension estimation using the box-counting method
Author: Aleksandra Pierzga
Context: Master's thesis project (2025)
Libraries used: porespy, skimage.io, skimage.color, numpy, os
"""

# Import libraries
import porespy as ps
from skimage.io import imread
from skimage.color import rgb2gray
import numpy as np
import os


def wymiar_fraktalny(sciezka, prog):
    '''
    This function estimates the fractal dimension of an image using the box-counting method.

    Parameters:
    sciezka - path to the image file
    prog - image binarization threshold value

    Returns:
    fd - fractal dimension estimated from the slope of the log–log relationship
    '''
    img = imread(sciezka)  # Load image
    if img.shape[-1] == 4:
        img = img[:, :, :3]  # Remove alpha (transparency) channel
    bw = rgb2gray(img)  # Convert image to grayscale
    binarnosc = bw < prog  # Pixels darker than the threshold are True, others are False
    wynik = ps.metrics.boxcount(binarnosc)  # Compute fractal dimension
    fd = np.mean(wynik['slope'])  # Compute the mean log–log slope
    return fd


folder = input("Enter the path to the folder with images: ")  # Read directory path
prog = float(input("Enter the binarization threshold in %: ")) / 100  # Enter threshold in percent

# Process all files in the directory
for plik in os.listdir(folder):
    if plik.lower().endswith(('.png', '.jpg', '.jpeg')):  # Filter for image files
        sciezka = os.path.join(folder, plik)
        fd = wymiar_fraktalny(sciezka, prog)
        print(plik, ":", fd)  # Print filename and result
