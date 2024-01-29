#@ File (label = "Project directory", style = "directory") proj
#@ String (label = "Previous step", choices = {"01b", "02"}, style = "radioButtonHorizontal") previous_step

// Script info
pipeline = "cell_segmentation_cellpose_measurement";
version = "v3.3";

// Input/Output
input = proj + File.separator + "01_auto_segmentation";
if (!File.exists(input)) exit("Input directory 01_auto_segmentation does not exist.");
output = input;
if (previous_step == "01b") {
	suffix = "_masks_cellpose.png";
	out_suffix = "_ROIs_cellpose.zip";
}
else if (previous_step == "02") {
	suffix = "_masks_auto_manual.png";
	out_suffix = "_ROIs_auto_manual.zip";
}

// Parameters
max_cellpose_artifact_roi_area = 10;  // px

// Main
close_everything();
processFolder(input);
close_everything();

function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		file = list[i];
		if (endsWith(file, suffix)) {
			// Extract sample name
			sample = rstrip(file, suffix);
			// Analyze image if this sample has not been completely analyzed (i.e., its output file does not exist)
			if (!isin(sample + out_suffix, getFileList(output))) {
				processFile(input, output, file);
			}
		}
	}
}

function processFile(input, output, file) {
	
	// Open downsized image
	// Why: looping through ROIs is somehow much faster on this tif image than on the png image
	open(input + File.separator + sample + "_downsize.tif");
	
	// Open Cellpose label image and convert to ROIs with filtering
	// Cellpose can create tiny disconnected objects: same pixel value as a closeby cell, often 1px big and sometimes 2px (maybe more), can be one or more (I have seen 1-3) per real cell
	// Each cell should be a connected component, thus convert each connected component into an ROI (instead of all pixels with the same pixel value into an ROI which could be a composite ROI) then filter out tiny ROIs
	open(input + File.separator + file);
	run("Label image to ROIs");  // very fast but cannot be used in batch mode 'hide'
	close();
	rn = roiManager("count");
	for (r = rn-1; r >= 0; r--) {
		roiManager("select", r);
		if (getValue("Area raw") < max_cellpose_artifact_roi_area) roiManager("delete");
	}
	clear_selection();
	
	// Rename ROIs
	rename_rois();
	// Save ROIs
	roiManager("Deselect");
	roiManager("Save", output + File.separator + sample + out_suffix);

}
