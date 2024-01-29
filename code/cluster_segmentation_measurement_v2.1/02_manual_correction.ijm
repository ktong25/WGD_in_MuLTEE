#@ File (label = "Project directory", style = "directory") proj
#@ String (label = "Input image file suffix", value = "_downsize.tif") suffix
#@ Boolean (label = "Do manual correction of auto segmentation results?", value = true) manual_correction
#@ String (label = "Grid interval (um) and color (e.g., '2000 Cyan') (if blank then not show grid)", value = "") grid_pars  // if "" then do not show grid to aid manual correction
#@ Float (label = "Min cluster area (um^2) for filtering after manual correction", value = 40) min_cluster_area
#@ Boolean (label = "Use batch mode 'hide' (if skip manual correction)", value = false) batch_mode_hide

// Check input parameter values
if (manual_correction && batch_mode_hide) exit("Cannot use batch mode 'hide' if do manual correction.");

// Script info
pipeline = "cluster_segmentation_measurement";
version = "v2.1";
step = "02";
task = "manual_correction";
script = String.join(newArray(pipeline, version, step + "_" + task), " ");

// Input/Output
input = proj + File.separator + "01_auto_segmentation";
output = proj + File.separator + "02_manual_correction";
if (!File.exists(input)) exit("Input directory 01_auto_segmentation does not exist.");
if (!File.exists(output)) File.makeDirectory(output);

// Print and save log
print_and_save_log("--------------------"); 
print_and_save_log("");
print_and_save_log("Running " + script + " ...");
print_and_save_log("ImageJ version: " + getVersion());
print_and_save_log("OS: " + getInfo("os.name"));
print_and_save_log("Project directory: " + proj);
print_and_save_log("Input directory: " + input);
print_and_save_log("Output directory: " + output);
print_and_save_log("Input image file suffix: " + suffix);
print_and_save_log("Do manual correction of auto segmentation results?: " + tf(manual_correction));
print_and_save_log("Grid interval (um) and color (if blank then not show grid): " + grid_pars);
print_and_save_log("Min cluster area (um^2) for filtering after manual correction: " + min_cluster_area);
print_and_save_log("Use batch mode 'hide': " + tf(batch_mode_hide));
print_and_save_log("");

// Process input parameter values
grid_pars = split(grid_pars, " ");

// Initialization
run("Colors...", "foreground=white background=black selection=red");
run("Options...", "iterations=1 count=1 black pad");  // binary operations
if (batch_mode_hide) setBatchMode("hide");  // set batch mode

// Main
close_everything();
if (manual_correction) {
	waitForUser("Make sure all adjustable parameter values in manual correction macros are appropriate.");
}
processFolder(input);
print_and_save_log("--------------------"); 
close_everything();

function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		file = list[i];
		if (endsWith(file, suffix)) {
			// Extract sample name
			sample = rstrip(file, suffix);
			// Specify input and output ROIs zip file names
			in_ROIs_zip = sample + "_ROIs_auto.zip";
			out_ROIs_zip = sample + "_ROIs_auto_manual.zip";
			out_ROIs_eoe_zip = sample + "_ROIs_auto_manual_eoe.zip";
			// Analyze image if this sample has not been completely analyzed (i.e., its final output file does not exist)
			if (!isin(out_ROIs_eoe_zip, getFileList(output))) {
				processFile(input, output, file);
			}
		}
	}
}

function processFile(input, output, file) {

	// Print sample info
	print_and_save_log("Processing sample " + sample + " ...");
	
	// Open image (downsized) for manual correction
	// This is necessary even if skip manual correction, in order to exclude cluster ROIs on image edges later
	open(input + File.separator + file);
	// Check if unit is micron (important for area thresholds, etc.)
	getPixelSize(unit, pixelWidth, pixelHeight);
	if (unit != "microns") exit("Unit is not microns");
	if (pixelWidth != pixelHeight) exit("pixelWidth != pixelHeight");

	// Load well center ROI if present
	if (isin(sample + "_wellcenter.roi", getFileList(input))) {
		roiManager("Open", input + File.separator + sample + "_wellcenter.roi");
		if (isin(sample + "_wellcenter_midline.roi", getFileList(input))) {  // if only include clusters touching left side of well center border
			roiManager("Open", input + File.separator + sample + "_wellcenter_midline.roi");
		}
		roiManager("Deselect");
		roiManager("Set Color", "yellow");
		run("From ROI Manager");  // add ROI to overlay
		roiManager("Deselect"); roiManager("Delete");
	}

	// Load auto segmentation ROIs for manual correction
	// Check and ask if continue from last time
	if (isin(out_ROIs_zip, getFileList(output))) {
		if (getBoolean("A pre-existing ROIs_auto_manual.zip is detected\n" +
					   "Continue from last time by loading ROIs_auto_manual.zip: click Yes\n" +
					   "Start from scratch by loading ROIs_auto.zip: click No\n" +
					   "Abort macro: click Cancel")) {
			print_and_save_log("Continue from last time ...");
			roiManager("Open", output + File.separator + out_ROIs_zip);
		}
		else {
			print_and_save_log("Start from scratch rather than continue from last time ...");
			roiManager("Open", input + File.separator + in_ROIs_zip);
		}
	}
	else roiManager("Open", input + File.separator + in_ROIs_zip);
	// Report cluster number
	print_and_save_log("Cluster number at input: " + roiManager("count"));
	
	// Wait for user to see the image and ROIs if no manual correction and not batch mode hide
	if ((!manual_correction) && (!batch_mode_hide)) {
		roiManager("show all without labels"); 
		waitForUser("Check and click OK to move to the next image");
	}
	
	// Manual correction
	if (manual_correction) {
		// Show grid as overlay to aid manual correction
		if (grid_pars.length > 0) run("Grid...", "grid=Lines area=" + pow(parseFloat(grid_pars[0]),2) + " color=" + grid_pars[1] + " center");
		// Manual correction
		roiManager("Show All without labels");
		correct_complete = false;
		while (!correct_complete) {
			waitForUser("Manual correction: \n" +
						" \n" +
						"Keyboard shortcuts to ImageJ commands: \n" +
						"- zoom tool [b], scrolling tool [n] \n" +
						"- rectangle tool [f], freehand selection tool [g], freeline selection tool [h] \n" +
						"- ROI manager: show all with labels [j], show none [k], show all without labels [l] \n" +
						" \n" +
						"Operations: \n" +
						"- Select an ROI: [j] + [g/h] + single click ROI number, or [j/l] + [g] + double click inside ROI \n" +
						"- Draw an undetected object: [g] + draw + [t] \n" +
						"- Delete an incorrect object: [g] + select + [x] (in practice: you can hold [x] and select ROIs one by one) \n" +
						"- Correct a detected object: [g] + draw + [t] + select + add [o] or remove [p] a region \n" +
						"- Merge two objects by simply combining them: [g] + select A + [e] + select B + [o] \n" +
						"- Merge two objects by filling 1px gap between them: [g] + select A + [e] + select B + [y] \n" +
						"- Split touching objects: [g/h] + draw + [t] + select + [u] \n" +
						"- Delete small objects in the whole image: [q] \n" +
						"- Delete all objects (entirely) within a selection: [g] + draw + [w] \n" +
						"- Auto segment cluster(s) in a selection (based on auto thresholding): [g] + draw + [r] \n" +
						"   - note: in the selection, include some background area, and shading should not be too variable \n" +
						"- Auto split touching clusters (based on marker-controlled watershed): [g] + select + [i] \n" +
						" \n" +
						"Recommendations for macroscopic strains: \n" +
						"- pipetting/sampling macroscopic clusters (> ~700um in diameter) is stochastic \n" +
						"- aim to detect (well-center-area/well-area)x100% (e.g. 70%) of macroscopic clusters in the whole well \n" +
						"- if fewer than 5 macroscopic clusters in the whole well, then you may detect all of them \n" +
						" \n" +
						"Faster workflow for some macroscopic strains with far more small debris than small clusters: \n" +
						"1. Preliminarily assign all ROIs as cluster (red) or debris (white) based on an area threshold [Q] \n" +
						"2. Manually go through all white ROIs and assign/rescue clusters (false negatives) as red [1] \n" +
						"    - If you make a mistake, you can assign ROIs back as white [2] \n" +
						"3. Delete all white ROIs [Z] (step 1-3 essentially remove all small debris in a fast way) \n" +
						"4. Perform manual correction (for all red ROIs) as you normally do \n" +
						"[Optional] For special purposes like generating cluster-debris annotations for machine learning: \n" +
						"1-2. Follow step 1-2 above \n" +
						"3. Manually go through all red ROIs and assign debris (false positives) as white [2] \n" +
						"4. Save all ROIs as e.g. *_ROIs_auto_manual_cluster_debris.zip (which saves ROI colors) in a folder \n" +
						"5-6. Follow step 3-4 above \n" +
						" \n" +
						"Click OK when done.\n" +
						"If somehow the macro aborts when correcting an image, manually save current ROIs as <sample>_ROIs_auto_manual.zip");
			correct_complete = getBoolean("Have you completed manual correction?");  // in case of error during manual correction, script proceeds to execute this line, thus user can still go back to continue manual correction and save all previous efforts
		}
		// Report cluster number
		print_and_save_log("Cluster number after manual correction: " + roiManager("count"));
		// Remove manual correction artifacts
		/// Mistakenly use freehand line tool to draw an ROI
		rn = roiManager("count");
		for (r = rn-1; r >= 0; r--) {
		    roiManager("select", r);
			if (getValue("Length") > 0) {
				run("To Selection");  // zoom in to selection
				run("Out [-]"); run("Out [-]"); run("Out [-]"); run("Out [-]"); run("Out [-]");
				roiManager("show all with labels");
				setTool("freehand");
				waitForUser("This is a freeline ROI.\n" +
							"Delete it, re-draw by freehand selection tool and add to ROI manager.\n" + 
							"Click OK when done.");
			}
		}
		clear_selection();
		/// Remove tiny objects (e.g., due to Split command in "Split ROI" macro)
		rn = roiManager("count");
		for (r = rn-1; r >= 0; r--) {
		    roiManager("select", r);
			if (getValue("Area") < min_cluster_area) roiManager("delete");
		}
		clear_selection();
		// Report cluster number
		print_and_save_log("Cluster number after removing manual correction artifacts: " + roiManager("count"));
		// Rename ROIs
		rename_rois();
	}
	
	// Remove overlay if there is any
	run("Remove Overlay");
	
	// Save ROIs
	roiManager("Deselect");
	roiManager("Save", output + File.separator + out_ROIs_zip);

	// Remove clusters on image edges
	rn = roiManager("count");
	for (r = rn-1; r >= 0; r--) {
		roiManager("select", r);
		getSelectionBounds(x, y, w, h);
		if ( (x<=0) || (y<=0) || (x+w>=getWidth()) || (y+h>=getHeight()) ) roiManager("delete");
	}
	clear_selection();
	// Report cluster number
	print_and_save_log("Cluster number after excluding clusters on edges: " + roiManager("count"));
	// Rename ROIs
	rename_rois();
	// Save ROIs
	roiManager("Deselect");
	roiManager("Save", output + File.separator + out_ROIs_eoe_zip);
	
	// Close
	// Close all images
	run("Close All"); 
	// Close ROI manager
	roiManager("deselect"); roiManager("delete"); if (!batch_mode_hide) {selectWindow("ROI Manager"); run("Close");}
	
	// Print complete
	print_and_save_log("Complete processing.");
	print_and_save_log("");	

}
