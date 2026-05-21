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

from tissue_analyzer import detect_nuclei, segment_marker

reader = BioImage("../data/cropped_for_testing/10389_Calnexin_001.czi")

nuclear_image = reader.data[0, 0, ...].squeeze()
marker_image = reader.data[0, 1, ...].squeeze()

tissue_mask = np.ones_like(marker_image, dtype=np.bool)

# plt.imshow(marker_image)
# plt.show()

marker_mask = segment_marker(marker_image, tissue_mask, threshold_mult=7)

# ov = sk.segmentation.mark_boundaries(sk.exposure.rescale_intensity(marker_image, out_range=(0.0, 1.0)), marker_mask)

# plt.imshow(ov)
# plt.show()

# plt.imsave("../processed/2026-05-21 Dev/10389_Calnexin_001_labeled.png", ov)

nuclear_mask = detect_nuclei(nuclear_image)
