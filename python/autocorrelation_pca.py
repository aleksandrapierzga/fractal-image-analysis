"""
Autocorrelation and PCA analysis of lichen images
Author: Aleksandra Pierzga
Context: Master's thesis project (2025)
Libraries used: os, numpy, matplotlib.pyplot, skimage.io, skimage.color, scipy.fft, sklearn.decomposition
"""

# Import libraries
import os
import numpy as np
import matplotlib.pyplot as plt
from skimage.io import imread
from skimage.color import rgb2gray
from scipy.fft import fft2, ifft2
from sklearn.decomposition import PCA


def profil_radialny(dane, srodek):
    '''
    This function computes the mean pixel value as a function of distance from the image center.

    Parameters:
    dane - 2D array (matrix)
    srodek - image center coordinates (y, x)

    Returns:
    Radial profile - a vector of mean pixel values for successive distances from the center.
    '''

    y, x = np.indices(dane.shape)  # Create a coordinate grid with the same shape as the autocorrelation map
    r = np.sqrt((x - srodek[1])**2 + (y - srodek[0])**2).astype(int)  # Distance from the center
    suma_wartosci = np.bincount(r.ravel(), dane.ravel())  # Sum of pixel values per radius
    liczba_punktow = np.bincount(r.ravel())
    srednie = suma_wartosci / liczba_punktow  # Mean pixel value for each radius
    return srednie


def autokorelacja_radialna(sciezka):
    '''
    This function computes the radial autocorrelation profile of an image using the Fast Fourier Transform (FFT)
    and displays the 2D autocorrelation map.

    Parameters:
    sciezka - path to the image

    Returns:
    A 1D correlation profile as a function of radius.
    '''
    img = imread(sciezka)  # Load image from file
    if img.shape[-1] == 4:  # Remove alpha (transparency) channel
        img = img[:, :, :3]
    bw = rgb2gray(img)  # Convert to grayscale
    f_bw = fft2(bw)
    auto = np.real(ifft2(f_bw * np.conj(f_bw)))  # Compute autocorrelation on a torus using the Fourier transform
    auto = np.fft.fftshift(auto)

    centrum = (auto.shape[0] // 2, auto.shape[1] // 2)
    profil = profil_radialny(auto, centrum)

    # 2D autocorrelation plots
    plt.figure(figsize=(6, 6))
    plt.imshow(auto, cmap='inferno')
    plt.colorbar()
    plt.title("Image autocorrelation")
    plt.axis('off')
    plt.tight_layout()
    plt.show()

    return profil


def analiza_pca_autokorelacji(foldery):
    '''
    This function performs PCA on previously computed radial autocorrelation profiles
    and displays the results as a plot.

    Parameters:
    foldery - a dictionary with species names and corresponding paths

    Returns:
    Displays a PCA plot.
    '''

    # Process images from the folders
    dane = []
    for gatunek, folder in foldery.items():
        for plik in os.listdir(folder):
            if plik.lower().endswith(('.png', '.jpg', '.jpeg')):  # Process image files only
                sciezka = os.path.join(folder, plik)
                profil = autokorelacja_radialna(sciezka)
                dane.append({
                    "gatunek": gatunek,
                    "plik": plik,
                    "profil": profil
                })

    # Determine the minimum profile length
    dlugosci = []
    for obiekt in dane:
        dlugosci.append(len(obiekt["profil"]))
    dlugosc = min(dlugosci)

    # Trim profiles to the same length and build matrix X
    X = []
    for obiekt in dane:
        X.append(obiekt["profil"][:dlugosc])
    X = np.array(X)

    # Build the list of species labels
    gatunki = []
    for obiekt in dane:
        gatunki.append(obiekt["gatunek"])

    # PCA
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X)

    # PCA plot
    plt.figure(figsize=(8, 6))
    gatunki_unikalne = sorted(set(gatunki))
    kolory = ['red', 'green', 'blue']

    for i, gatunek in enumerate(gatunki_unikalne):
        indeksy = [j for j, g in enumerate(gatunki) if g == gatunek]
        plt.scatter(X_pca[indeksy, 0], X_pca[indeksy, 1], label=gatunek, color=kolory[i % len(kolory)], alpha=0.7)

    plt.title("PCA of autocorrelation profiles")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend()
    plt.grid(True)
    plt.show()


# Paths to folders
folders = {}
gatunki = ["Cladonia rangiferina", "Hypogymnia physodes", "Xanthoria parietina"]

for gatunek in gatunki:
    sciezka = input("Enter the path to the folder with images for '{}': ".format(gatunek))
    folders[gatunek] = sciezka

# Function call
analiza_pca_autokorelacji(folders)
