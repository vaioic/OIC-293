import csv
import re
from pathlib import Path

import numpy as np
import skimage as sk
from bioio import BioImage
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.widgets import RectangleSelector
from scipy import ndimage


def process_directory(input_dir, output_dir):

    # Validate the inputs
    if isinstance(input_dir, str):
        input_dir = Path(input_dir)
    elif isinstance(input_dir, Path):
        pass
    else:
        raise TypeError(f"Expected input_dir to be a str or Path. Instead it is a {type(input_dir)}.")
    
    if not input_dir.exists():
        input_dir.mkdir(parents=True)
    elif input_dir.is_file():
        raise TypeError(f"Expected input_dir to be a directory. Instead it appears to be a file: {input_dir}.")
    
    if isinstance(output_dir, str):
        output_dir = Path(output_dir)
    elif isinstance(output_dir, Path):
        pass
    else:
        raise TypeError(f"Expected output_dir to be a str or Path. Instead it is a {type(output_dir)}.")
    
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
    elif output_dir.is_file():
        raise TypeError(f"Expected output_dir to be a directory. Instead it appears to be a file: {output_dir}.")
    

    with open(output_dir / "measurements.csv", "w", newline="") as csvfile:

        # Write CSV headers
        csvwriter = csv.writer(csvfile, delimiter=",")
        csvwriter.writerow(["File", "Marker Area (px)", "Tissue Area (px)", "Ratio"])                          
    
        file_list = get_filtered_file_list(input_dir)

        for f in file_list:

            num_pixels_marker, num_pixels_tissue, ratio = process_image(f, output_dir)

            csvwriter.writerow([f.name, num_pixels_marker, num_pixels_tissue, ratio])    

def get_filtered_file_list(input_dir):

    file_list = input_dir.glob("*.czi")

    pattern = re.compile(r"^([0-9]+)_([\w]+)(-rescan)?\.*", re.IGNORECASE)

    filtered_files = {}

    for f in file_list:
        match = pattern.match(f.name)

        sample_id, marker, rescan_flag = match.groups()

        # To help keep files unique, we create a dict
        key = (sample_id, marker.lower())

        if (key not in filtered_files) or (rescan_flag is not None):
            filtered_files[key] = f

    # Extract the final list of Path objects
    final_file_list = list(filtered_files.values())
    
    return final_file_list

def process_image(image_path, output_path, config_path=None, downsample=0.25):

    # Validate the inputs
    if isinstance(image_path, str):
        image_path = Path(image_path)
    elif isinstance(image_path, Path):
        pass
    else:
        raise TypeError(f"Expected image_path to be a str or Path. Instead it is a {type(image_path)}.")
    
    if isinstance(output_path, str):
        output_path = Path(output_path)
    elif isinstance(output_path, Path):
        pass
    else:
        raise TypeError(f"Expected output_path to be a str or Path. Instead it is a {type(output_path)}.")
    
    if not output_path.exists():
        output_path.mkdir(parents=True)
    elif output_path.is_file():
        raise TypeError(f"Expected output_path to be a directory. Instead it appears to be a file: {output_path}.")
    
    if config_path:
        if isinstance(config_path, str):
            config_path = Path(config_path)
        elif isinstance(config_path, Path):
            pass
        else:
            raise TypeError(f"Expected config_path to be a str or Path. Instead it is a {type(config_path)}.")
        
        if not config_path.is_fileI():
            raise TypeError(f"Config file was not found at path: {config_path}.")
        
    # Read in the image
    reader = BioImage(image_path)

    # Define the different channels
    img_DAPI = reader.data[0, 0, ...].squeeze()
    img_marker = reader.data[0, 1, ...].squeeze()

    if not config_path:

        # Get the ROI
        rois = get_ROI(img_marker)

    # Downsample the image (if requested)
    if (not downsample is None) and (downsample > 0):
        img_marker = sk.transform.rescale(img_marker, downsample)
        img_DAPI = sk.transform.rescale(img_DAPI, downsample)

    # Segment the tissue
    tissue_mask = segment_tissue(img_DAPI)

    # ov = sk.segmentation.mark_boundaries(img_DAPI, tissue_mask, mode="thick")

    # plt.imshow(ov)
    # plt.show()
    # exit()

    # Segment the marker
    marker_mask = segment_marker(img_marker, tissue_mask)

    # Generate the output image
    output_img = generate_output_image(img_marker, tissue_mask, marker_mask, downsample=0.25)

    # Save the output image
    fn = image_path.stem
    sk.io.imsave(output_path / (fn + ".png"), output_img)

    # Measure the number of pixels (equiv. to area) of the tissue and marker
    num_pixels_tissue = np.count_nonzero(tissue_mask)
    num_pixels_marker = np.count_nonzero(marker_mask)
    ratio = num_pixels_marker / num_pixels_tissue

    return (num_pixels_marker, num_pixels_tissue, ratio)

def get_ROI(image):

    all_rois = []
    current_coords = None  # Holds the unsaved box coordinates

    def onselect(eclick, erelease):
        """Updates the temporary coordinates whenever a box is drawn/resized."""
        global current_coords
        xmin, xmax = int(min(eclick.xdata, erelease.xdata)), int(max(eclick.xdata, erelease.xdata))
        ymin, ymax = int(min(eclick.ydata, erelease.ydata)), int(max(eclick.ydata, erelease.ydata))
        current_coords = (xmin, xmax, ymin, ymax)

    def on_key(event):
        """Listens for keyboard inputs."""
        global current_coords
        
        if event.key == 'enter':
            if current_coords is not None:
                xmin, xmax, ymin, ymax = current_coords
                
                # 1. Save the coordinates permanently
                roi = {"xmin": xmin, "xmax": xmax, "ymin": ymin, "ymax": ymax}
                all_rois.append(roi)
                print(f"Saved ROI #{len(all_rois)}: {roi}")
                
                # 2. Draw a permanent visual patch on the image
                width = xmax - xmin
                height = ymax - ymin
                rect = Rectangle((xmin, ymin), width, height, edgecolor='red', facecolor='none', linewidth=2)
                ax.add_patch(rect)
                
                # 3. Refresh the plot to show the new permanent patch
                fig.canvas.draw()
                
                # Reset temporary storage so we don't duplicate on double-enter
                current_coords = None 
            else:
                print("No new ROI drawn to save!")

    # Downsize the image for easier viewing
    image_ds = sk.transform.rescale(image, 0.25)

    fig, ax = plt.subplots()
    ax.imshow(image_ds, cmap="gray")

    fig.canvas.mpl_connect('key_press_event', on_key)

    # Enable the selector
    rs = RectangleSelector(ax, onselect, useblit=True,
                           button=[1], interactive=True)

    plt.show()
    print(f"\nFinal saved ROIs count: {len(all_rois)}")
    # with open("multiple_rois.txt", "w") as f:
    #     for i, roi in enumerate(all_rois):
    #         f.write(f"ROI_{i}:{roi['xmin']},{roi['xmax']},{roi['ymin']},{roi['ymax']}\n")

    return all_rois 


def segment_tissue(image):

    #Filter the image
    image = sk.filters.gaussian(image, sigma=3)

    # Normalize the intensity
    image = sk.exposure.rescale_intensity(image, out_range=(0.0, 1.0))

    thresh = sk.filters.threshold_otsu(image)

    mask = image > (0.95 * thresh)

    mask = sk.morphology.remove_small_objects(mask, max_size=10000)
    mask = ndimage.binary_fill_holes(mask)
    #mask = sk.morphology.opening(mask, sk.morphology.disk(30))

    mask = sk.segmentation.clear_border(mask) 
    

    return mask

def segment_marker(image, tissue_mask, inset=100, threshold_mult=2):

    # Inset the tissue mask to avoid edge effects
    if inset is not None:
        tissue_mask = shrink_mask(tissue_mask, inset)

    # Normalize the intensity
    image = sk.exposure.rescale_intensity(image, out_range=(0.0, 1.0))

    #thresh = sk.filters.threshold_otsu(image[tissue_mask])
    #mask = image > (0.96 * thresh)

    #thresh = np.mean(image) + (1.5 * np.std(image))

    median_intensity = np.median(image)
    diff_from_median = np.abs(image - median_intensity)
    MAD = np.median(diff_from_median)

    scaled_MAD = 1.4826 * MAD

    thresh = median_intensity + (scaled_MAD)
    mask = image > thresh

    # plt.imshow(mask)
    # plt.show()
    # exit()

    # mask = sk.morphology.remove_small_objects(mask, max_size=1000)
    #mask = ndimage.binary_fill_holes(mask)

    mask[~tissue_mask] = False

    return mask

def shrink_mask(mask, shrink_by):

    distance = ndimage.distance_transform_edt(mask)
    mask[distance < shrink_by] = False

    return mask

def generate_output_image(image, tissue_mask, marker_mask, downsample=0.25):

    if not downsample is None:
        image = sk.transform.rescale(image, downsample)
        tissue_mask = sk.transform.rescale(tissue_mask, downsample)
        marker_mask = sk.transform.rescale(marker_mask, downsample)

    #image_norm = sk.exposure.equalize_hist(image)
    p_low, p_high = np.percentile(image, (5, 95))
    image_norm = sk.exposure.rescale_intensity(image, in_range=(p_low, p_high), out_range=(0.0, 1.0))

    # Insert the tissue outline
    output_img = sk.segmentation.mark_boundaries(image_norm, tissue_mask, mode="thick", color=(0, 1, 0))

    # Insert the identified marker outline
    output_img[..., 0] = 0.5 * output_img[..., 0] + 0.5 * marker_mask
    output_img[..., 1] = 0.5 * output_img[..., 1]
    output_img[..., 2] = 0.5 * output_img[..., 2] + 0.5 * marker_mask

    #output_img = sk.segmentation.mark_boundaries(output_img, marker_mask, mode="thick", color=(1, 0, 1))

    output_img = sk.util.img_as_ubyte(output_img)

    return output_img

def main():
    # This is primarily for testing

    # process_image("../data/10389_Plin2-rescan.czi", "../processed/2026-05-18 Dev")

    process_image("../data/10390_Plin2.czi", "../processed/2026-05-18 Dev")
    # process_directory("../data", "../processed/2026-05-18/")
    
    # process_directory("../data", "../processed/2026-05-15 Dev/")

if __name__ == "__main__":
    main()
