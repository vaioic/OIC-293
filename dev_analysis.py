import numpy as np
import skimage as sk
from bioio import BioImage
from matplotlib import pyplot as plt
from scipy import ndimage

reader = BioImage("../data/10389_Plin2-rescan.czi")

img = reader.data

# Channels are DAPI, and AF-555

# #T, C, Z, Y, X
# print(img.shape)

# plt.imshow((img[0, 1, ...]).squeeze())
# plt.show()

# Segment the organelle
img_ch2 = (img[0, 1, ...]).squeeze()

# Downsample the test image
img_ch2_ds = sk.transform.rescale(img_ch2, 0.25)

# Normalize the intensity
img_ch2_ds_norm = sk.exposure.rescale_intensity(img_ch2_ds, out_range=(0.0, 1.0))

#Filter the image
img_ch2_ds_norm = sk.filters.gaussian(img_ch2_ds_norm, sigma=1.5)

# First, find the tissue region
mask_tissue = img_ch2_ds_norm > 0.010

#mask_tissue = sk.morphology.opening(mask_tissue, sk.morphology.disk(30))
print(f"removing small objects")
mask_tissue = sk.morphology.remove_small_objects(mask_tissue, max_size=10000)

print(f"Filling holes")
mask_tissue = ndimage.binary_fill_holes(mask_tissue)

overlay_tissue = sk.segmentation.mark_boundaries(img_ch2_ds_norm, mask_tissue, color=(0, 1, 0))

# plt.imshow(overlay_tissue)
# plt.show()

# Now threshold the rest
thresh = sk.filters.threshold_otsu(img_ch2_ds_norm[mask_tissue])
mask = img_ch2_ds_norm > 0.95 * thresh

mask = sk.morphology.remove_small_objects(mask, max_size=1000)
mask = ndimage.binary_fill_holes(mask)

# Normalize the intensity
img_ch2_ds_output = sk.exposure.rescale_intensity(img_ch2_ds, out_range=(0.0, 1.0))


overlay = sk.segmentation.mark_boundaries(img_ch2_ds_output, mask, mode="thick", color=(1, 1, 0))

# plt.imshow(overlay)
# plt.show()

# Save images
sk.io.imsave("test_tissue.png", sk.util.img_as_ubyte(overlay_tissue))
sk.io.imsave("marked_region.png", sk.util.img_as_ubyte(overlay))

nnz_tissue = np.count_nonzero(mask_tissue)
nnz_protein = np.count_nonzero(mask)

# Return the area ratioi
print(f"Ratio of areas = {(nnz_protein/nnz_tissue) * 100}")
