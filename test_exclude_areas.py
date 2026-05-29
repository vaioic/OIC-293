#%% Load imports and image data
import csv
import json
import re
from datetime import datetime
from pathlib import Path

import numpy as np
import skimage as sk
from bioio import BioImage
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.widgets import RectangleSelector
from scipy import ndimage

import tissue_analyzer

image_path = "../data/10402_Plin2.czi"

# Read in the image
reader = BioImage(image_path)

# Define the different channels
img_DAPI = reader.data[0, 0, ...].squeeze()
img_marker = reader.data[0, 1, ...].squeeze()

print(img_marker.shape)

# Rescale image intensity for downstream processing
img_DAPI_norm = sk.exposure.rescale_intensity(img_DAPI, out_range=(0.0, 1.0))
img_marker_norm = sk.exposure.rescale_intensity(img_marker, out_range=(0.0, 1.0))

#%% Test tissue segmentation

# # Test shade correction
# img_crop = img_marker_norm[19429:24521, 30206:35669]

# img_marker_blurred = sk.filters.gaussian(img_marker_norm, sigma=30)
# img_marker_corrected = img_marker_blurred - img_marker_blurred

# plt.imshow(sk.exposure.equalize_hist(img_marker_corrected))
# plt.show()
# exit()

tissue_threshold = tissue_analyzer.get_threshold(img_marker_norm)

# Downscale the tissue image for faster processing
img_ds = img_marker_norm[::8, ::8]

tissue_mask = img_ds > tissue_threshold

tissue_mask = sk.morphology.remove_small_holes(tissue_mask, max_size=1000)
tissue_mask = sk.morphology.remove_small_objects(tissue_mask, max_size=8000)

# Shrink the mask
#$tissue_mask = sk.morphology.isotropic_erosion(tissue_mask, 300)

tissue_mask_solid = ndimage.binary_fill_holes(tissue_mask)
tissue_mask_final = sk.morphology.isotropic_erosion(tissue_mask_solid, 200)

tissue_mask_final[~tissue_mask] = False

#%% Plotting

#Overlay mask
alpha = 0.2

img_overlay = np.zeros((img_ds.shape[0], img_ds.shape[1], 3))

img_ds_brighter = sk.exposure.equalize_hist(img_ds)

img_overlay[..., 0] = (1 - alpha) * img_ds_brighter
img_overlay[..., 1] = (1 - alpha) * img_ds_brighter + (alpha * tissue_mask_final)
img_overlay[..., 2] = (1 - alpha) * img_ds_brighter

plt.imshow(img_overlay)
plt.show()

# %% Try to identify other abnormalities to exclude



