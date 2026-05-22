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


def process_directory(input_dir, output_dir):
    """
    Process a directory of images

    Parameters
    ----------
    input_dir : str or Path
        Path to the directory containing images to analyze
    output_dir : str or Path
        Path to the directory where output data will be saved

    Raises
    ------
    TypeError
        The input_dir must be a str or Path
    TypeError
        The input_dir appears to be a file instead of a directory path
    TypeError
        The output_dir must be a str or Path
    TypeError
        The output_dir appears to be a file instead of a directory path    
    """

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
    """
    Get a filtered file list.

    The dataset sometimes includes duplicate files. This function filters the dataset to avoid processing duplicate samples. In particular, if it encounters a filename with the suffix "-rescan", it keep this version over the original. The function also handles duplicate files if the filenames are in a different case.

    Parameters
    ----------
    input_dir : Path
        Path to directory containing images to analyze

    Returns
    -------
    filtered_file_list : list of Path
        The final filtered list of files
    """

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
    filtered_file_list = list(filtered_files.values())
    
    return filtered_file_list

def process_image(image_path, output_path, config_path=None, downsample=None, tissue_threshold=None, tissue_type="auto"):
    """
    Process an image

    Parameters
    ----------
    image_path : str or Path
        Path to the image file
    output_path : str or Path
        Path to directory to save data
    config_path : str or Path, optional
        Path to configuration file, by default None
    downsample : float, optional
        If set, the images are downsampled by this value, by default None

    Raises
    ------
    TypeError
        _description_
    TypeError
        _description_
    TypeError
        _description_
    TypeError
        _description_
    TypeError
        _description_
    """

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
    
    roi_list = []  # Initialize an empty list of ROIs to process

    if config_path:
        if isinstance(config_path, str):
            config_path = Path(config_path)
        elif isinstance(config_path, Path):
            pass
        else:
            raise TypeError(f"Expected config_path to be a str or Path. Instead it is a {type(config_path)}.")
        
        if not config_path.is_file():
            raise TypeError(f"Config file was not found at path: {config_path}.")
        
        # Load the configuration file
        with open(config_path, "r") as cf:
            config = json.load(cf)

        roi_list = config["rois"]
        downsample = config["settings"]["downsample"]

        print("Loaded configuration file.")
        print(roi_list)

    # Read in the image
    reader = BioImage(image_path)

    # Define the different channels
    img_DAPI = reader.data[0, 0, ...].squeeze()
    img_marker = reader.data[0, 1, ...].squeeze()

    # # TEMP:
    # img_DAPI = img_DAPI[::8, ::8]
    # img_marker = img_marker[::8, ::8]

    # Rescale image intensity for downstream processing
    img_DAPI_norm = sk.exposure.rescale_intensity(img_DAPI, out_range=(0.0, 1.0))
    img_marker_norm = sk.exposure.rescale_intensity(img_marker, out_range=(0.0, 1.0))

    if not tissue_threshold:
        # If the tissue_threshold value was not previously in the configuration file, estimate a new threshold
        print(f"Estimating tissue threshold value.")
        tissue_threshold = get_threshold(img_marker_norm)

    # Look at tissue type to determine threshold factor
    if (not tissue_type) or (tissue_type.lower() == "auto"):

        print(f"Getting tissue type")

        # Get tissue type from filename
        fn = (image_path.stem).lower()

        if "calnexin" in fn:
            tissue_type = "er"
        elif "plin2" in fn:
            tissue_type = "lipid"
        elif "tom20" in fn:
            tissue_type = "mito"
        else:
            raise ValueError(f"Could not determine tissue type from filename. Please specify it directly.")
        
    print(f"Tissue type: {tissue_type}")
            
    match tissue_type.lower():

        case "er" | "lipid":
            marker_threshold_factor = 7
        
        case "mito":
            marker_threshold_factor = 12

        case _:
            raise ValueError(f"Could not determine marker threshold factor for unknown tissue type.")
   
    print(f"Threshold factor: {marker_threshold_factor}")
    # exit()

    if not roi_list:
        # If no ROIs were provided in the configuration file, prompt the user to select some

        # Get ROIs to process
        roi_list = get_ROI(img_marker)

    # Parse each ROI to measure marked area
    all_data = []
    for idx, roi in enumerate(roi_list):

        # Index the sub-images
        roi_DAPI = img_DAPI_norm[roi["ymin"]:roi["ymax"], roi["xmin"]:roi["xmax"]]
        roi_marker = img_marker_norm[roi["ymin"]:roi["ymax"], roi["xmin"]:roi["xmax"]]

        data, tissue_mask, marker_mask, nucl_detections, nucl_labels = process_ROI(roi_DAPI, roi_marker, threshold=tissue_threshold, downsample=downsample, marker_threshold_factor=marker_threshold_factor)

        # Add ROI metadata
        data["filename"] = str(image_path)
        data["ROI_xmin"] = roi["xmin"]
        data["ROI_xmax"] = roi["xmax"]
        data["ROI_ymin"] = roi["ymin"]
        data["ROI_ymax"] = roi["ymax"]

        all_data.append(data)
        
        # Generate and save output images
        rgb_img = np.zeros((roi_DAPI.shape[0], roi_DAPI.shape[1], 3))

        # Normalize the intensities
        roi_marker = sk.exposure.rescale_intensity(roi_marker, out_range=(0.0, 1.0))
        roi_DAPI = sk.exposure.rescale_intensity(roi_DAPI, out_range=(0.0, 1.0))

        rgb_img[..., 0] = roi_marker
        rgb_img[..., 1] = roi_marker
        rgb_img[..., 2] = roi_DAPI

        # plt.imshow(rgb_img)
        # plt.show()

        ov_mask = sk.segmentation.mark_boundaries(rgb_img, tissue_mask, 
                                                mode="thick", color=(0, 1, 0))
        ov_mask = sk.segmentation.mark_boundaries(ov_mask, marker_mask, 
                                                  mode="thick", 
                                                  color=(1, 0, 1))
        
        sk.io.imsave(output_path / ("roi" + f"{idx}" + "_masks.png"), sk.util.img_as_ubyte(ov_mask))

        # Mark the nuclei
        ov_nucl = sk.segmentation.mark_boundaries(
            roi_DAPI,
            nucl_labels)
        
        plt.imshow(ov_nucl)

        # Get positions of nuclei
        xx = []
        yy = []

        for d in nucl_detections:

            y, x, _ = d
            xx.append(x)
            yy.append(y)

        plt.scatter(xx, yy, 1)
        plt.savefig(output_path / ("roi" + f"{idx}" + "_nucl.png"))
                                                
        
    # # Save the outputs
    fn = image_path.stem  # Prefix for saved files

    # # Generate the output image
    # output_img = generate_output_image(img_marker, tissue_mask, marker_mask, downsample=0.25)

    # # Save the output image
    # sk.io.imsave(output_path / (fn + ".png"), output_img)

    # Save data
    all_keys = set().union(*(d.keys() for d in all_data))

    headers = ["filename", "ROI_xmin", "ROI_xmax", "ROI_ymin", "ROI_ymax"]

    # Automatically append the rest of the keys
    for key in sorted(all_keys):
        if key not in headers:
            headers.append(key)
    
    with open(output_path / (fn + "_data.csv"), mode="w", newline="", encoding="utf-8") as file:
            
        writer = csv.DictWriter(file, fieldnames=headers)
        
        writer.writeheader()  # Writes the first row (column names)
        writer.writerows(all_data)  # Writes all rows of data

    # Save processing configuration
    config_data = {
        "metadata": {
            "date_created": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "file": str(image_path),
            "config_version": "1.0"
        },
        "settings": {
            "downsample": downsample
        },
        "rois": roi_list
    }

    with open(output_path / (fn + "_config.json"), "w") as cf:
        json.dump(config_data, cf, indent=4)

    print(f"Config file saved to {str(output_path / (fn + '_config.json'))}")

    # # Measure the number of pixels (equiv. to area) of the tissue and marker
    # num_pixels_tissue = np.count_nonzero(tissue_mask)
    # num_pixels_marker = np.count_nonzero(marker_mask)
    # ratio = num_pixels_marker / num_pixels_tissue

    # return (num_pixels_marker, num_pixels_tissue, ratio)

def process_ROI(img_DAPI, img_marker, threshold=None, downsample=None, marker_threshold_factor=3):

    # Downsample the image (if requested)
    if (not downsample is None) and (downsample > 0):
        img_marker = sk.transform.rescale(img_marker, downsample)
        img_DAPI = sk.transform.rescale(img_DAPI, downsample)
        
    # Segment the tissue
    tissue_mask = segment_tissue(img_marker, threshold=threshold)

    if not np.any(tissue_mask):
        raise ValueError(f"No tissue section was found.")
    
    # Segment the marker
    marker_mask = segment_marker(img_marker, tissue_mask, threshold_mult=marker_threshold_factor)

    # Identify nuclei - Detections are likely more accurate for number of nuclei, as some do not show up in the labels
    nucl_detections, nucl_labels = detect_nuclei(img_DAPI)

    nucl_props = sk.measure.regionprops_table(
        nucl_labels, 
        properties=["area", "eccentricity"]
    )

    # Measure data
    data = {
        "num_nuclei": len(nucl_detections),
        "mean_nuclear_area": np.mean(nucl_props["area"]),
        "mean_nuclear_circularity": np.mean(nucl_props["eccentricity"]),
        "total_tissue_area": np.count_nonzero(tissue_mask),
        "total_marker_area": np.count_nonzero(marker_mask),
        "ratio_areas": np.count_nonzero(marker_mask) / np.count_nonzero(tissue_mask)
    }

    # ov = sk.segmentation.mark_boundaries(img_DAPI, nuclear_mask)

    # plt.imshow(ov)
    # plt.show()
    return data, tissue_mask, marker_mask, nucl_detections, nucl_labels

def detect_nuclei(image):

    image = sk.filters.gaussian(image, 2)

    threshold = sk.filters.threshold_otsu(image)

    mask = image > (threshold)
    mask = sk.morphology.remove_small_holes(mask, max_size=500)

    # Find centers for watershed
    blobs = sk.feature.blob_log(image, min_sigma=12, max_sigma=30, threshold=0.000005)

    # Filter the blobs by intensity to 

    # Make the markers
    markers = np.zeros_like(mask, dtype=np.int32)

    for idx, blob in enumerate(blobs):
        y, x, c = blob
        markers[int(y), int(x)] = idx + 1

    dd = ndimage.distance_transform_edt(mask)
    labels = sk.segmentation.watershed(-dd, markers, mask=mask, compactness=0.5)

    ov = sk.segmentation.mark_boundaries(sk.exposure.equalize_hist(image), labels)

    plt.imshow(ov)
    plt.show()

    return blobs, labels

def get_ROI(image, downsample_factor=8):

    all_rois = []
    current_coords = None  # Holds the unsaved box coordinates

    def onselect(eclick, erelease):
        """Updates the temporary coordinates whenever a box is drawn/resized."""
        nonlocal current_coords
        
        xmin, xmax = int(min(eclick.xdata, erelease.xdata)), int(max(eclick.xdata, erelease.xdata))
        ymin, ymax = int(min(eclick.ydata, erelease.ydata)), int(max(eclick.ydata, erelease.ydata))
        current_coords = (xmin, xmax, ymin, ymax)

    def on_key(event):
        """Listens for keyboard inputs."""
        nonlocal current_coords
                
        if event.key == 'enter':
            if current_coords is not None:
                xmin, xmax, ymin, ymax = current_coords
                
                # Save the coordinates as a dict
                roi = {"xmin": xmin, "xmax": xmax, "ymin": ymin, "ymax": ymax}
                all_rois.append(roi)
                print(f"ROI #{len(all_rois)}: {roi}")
                
                # Draw a rectangle on the image
                width = xmax - xmin
                height = ymax - ymin
                rect = Rectangle((xmin, ymin), width, height, edgecolor='green', facecolor='none', linewidth=1)
                ax.add_patch(rect)
                
                # Refresh the plot to show the new permanent patch
                fig.canvas.draw()
                
                # Reset temporary storage so we don't duplicate on double-enter
                current_coords = None 
            else:
                print("No new ROI drawn to save!")

    # Downsize the image for easier viewing
    if downsample_factor:
        image = image[::downsample_factor, ::downsample_factor]

    image = sk.exposure.rescale_intensity(image, out_range=(0.0, 1.0))

    fig, ax = plt.subplots(figsize=(12, 10))
    ax.imshow(image, cmap="gray")
    ax.set_title("Drag to resize and move the selection. Press enter to create an ROI. Close image when done.")

    fig.canvas.mpl_connect('key_press_event', on_key)

    # Enable the selector
    rs = RectangleSelector(ax, onselect, useblit=True,
                           button=[1], interactive=True)

    plt.show()

    # Rescale the ROIs by the downsample factor
    if downsample_factor:

        final_roi_list = [
            {
                "xmin": roi["xmin"] * downsample_factor,
                "xmax": roi["xmax"] * downsample_factor,
                "ymin": roi["ymin"] * downsample_factor,
                "ymax": roi["ymax"] * downsample_factor,
            }
            for roi in all_rois
        ]

    else:
        final_roi_list = all_rois

    return final_roi_list

def get_threshold(image, ds=16):

    thresh = sk.filters.threshold_otsu(image[::ds, ::ds])
   
    return thresh

def segment_tissue(image, threshold):

    #Filter the image
    #image = sk.filters.gaussian(image, sigma=3)


    if not threshold:
        thresh = get_threshold(image)

    mask = image > threshold

    mask = sk.morphology.remove_small_objects(mask, max_size=10000)
    #mask = ndimage.binary_fill_holes(mask)
    #mask = sk.morphology.opening(mask, sk.morphology.disk(30))

    #mask = sk.segmentation.clear_border(mask)     

    return mask

def segment_marker(image, tissue_mask, inset=None, threshold_mult=7):

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

    thresh = median_intensity + (threshold_mult * scaled_MAD)
    mask = image > thresh

    # plt.imshow(mask)
    # plt.show()
    # exit()

    mask = sk.morphology.remove_small_objects(mask, max_size=20)
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

    #process_image("../data/10390_Plin2.czi", "../processed/2026-05-18 Dev")
    process_image("../data/10390_Plin2.czi", "../processed/2026-05-22 Dev")
    # process_directory("../data", "../processed/2026-05-18/")    
    # process_directory("../data", "../processed/2026-05-15 Dev/")

if __name__ == "__main__":
    main()
