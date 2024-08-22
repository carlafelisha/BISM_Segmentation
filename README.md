# BISM_Segmentation.m

## Overview
This MATLAB script processes Traction Force Microscopy (TFM) data and corresponding spheroid image sequences to perform the following tasks:

1. **BISM Calculations**: Code referenced from_ Nier, V. et al., Inference of internal stress in a cell monolayer, Biophysical Journal, 110(7), 1625-1635, 2016_
2. **Spheroid Segmentation**: Identify the outline of the spheroid and create masks for different regions surrounding it (Regions A: 0-10μm, B: 10-20μm, C: -10-0μm, where 0 represents the edge of the spheroid).
3. **Stress Field Extraction**: Extract the stress field values in the defined regions around the spheroid.
4. **Averaging of Values**: Calculate the average tensile and compressive stress values for the different regions at each timeframe.

## Inputs
- **Spheroid Image Sequence**: Folder containing a sequence of JPEG images of the spheroid.
- **TFM Data**: Traction Force Microscopy data in PIV format.

## Outputs
- **Spheroid Mask Images**: Images displaying the spheroid mask to verify the segmentation process.
- **Stress Values Excel Sheet**: An Excel file containing the average tensile and compressive stress values for different regions at each timeframe. This file is saved in the directory where the MATLAB script is located.

## Workflow and Setup

### Step 1: Preparing the Environment
1. **Convert TIFF to JPG**: Convert your spheroid image sequence from TIFF to JPG format using ImageJ/Fiji. 
   - Open the TIFF file in ImageJ.
   - Adjust brightness/contrast: `Image > Adjust > Brightness/Contrast`, then select "Auto" and "Apply."
   - Save the sequence as JPG: `File > Save As > Image Sequence`, choose JPEG format, and specify the folder.

2. **Modify Directory Paths**:
   - Update the directory for the JPG images folder:
     - `[Line 7] dirPath`
   - Update the directory for the TFM FTTC/PIV files:
     - `[Line 18 and 137] ForceName`
     - `[Line 10] DirectoryName`
     - `[Line 11] Sessions`

3. **Set Coefficient Based on Magnification**:
   - `[Line 27/28] coeff = 0.325` for 40x magnification or `0.216` for 60x magnification.

4. **Adjust Timeframes**:
   - Modify `k0` and `kM` to set the start and end timeframes for processing:
     - `[Line 49] kM`

### Step 2: Image Segmentation
1. **Adjust File Naming**: Modify the filename based on the naming convention used during the JPG conversion:
   - `[Line 344] fileName = sprintf('%s_dic%03d.jpg', Session, k0);`

2. **Tweak Segmentation Parameters**:
   - Fine-tune the segmentation parameters, such as the "radius" for erosion and dilation, to optimize the detection of the spheroid. 
   - The erosion process removes small details, while dilation highlights important features.
   - Check the segmentation results for a few timeframes and adjust the erosion and dilation parameters as needed.
   - Read more about these operations here: [Types of Morphological Operations (MATLAB)](https://www.mathworks.com/help/images/morphological-dilation-and-erosion.html)
 

### Limitations
- The segmentation roughly selects the outline of the largest detectable object in the frame.
- It cannot automatically detect spheroid breakdown/spreading; data from those frames must be manually removed.
- The generated segmentation masks are resized to fit the resolution of the BISM output. Thus, there is loss of detail, potentially reducing the accuracy of the masks in capturing fine structural features.
- BISM stress field graphs should also be plotted to identify valid data points.
  - To generate Stress Field Images, lines 459-491 can be uncommented to output spheroid images overlaid with stress field heatmap
    
