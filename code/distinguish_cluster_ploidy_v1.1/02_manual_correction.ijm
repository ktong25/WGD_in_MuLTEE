#@ File (label = "Project directory", style = "directory") proj
#@ String (label = "Input image file suffix", value = "_downsize.tif") suffix
#@ Boolean (label = "Do manual correction of auto segmentation results?", value = true) manual_correction
#@ String (label = "Grid interval (um) and color (e.g., '2000 Cyan', '')", value = "") grid_pars  // if "" then do not show grid to aid manual correction
#@ Float (label = "Min ROI area (um^2) for filtering after manual correction", value = 0) min_roi_area
#@ Boolean (label = "Use batch mode 'hide' (if not do manual correction)", value = false) batch_mode_hide

// Check input parameter values
if (manual_correction & batch_mode_hide) exit("Cannot use batch mode 'hide' if do manual correction.");

// Script info
pipeline = "cell_segmentation_cellpose_measurement";
version = "v1.1";
step = "02_manual_correction";
script = String.join(newArray(pipeline, version, step), "_");

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
print_and_save_log("Project directory: " + proj);
print_and_save_log("Input directory: " + input);
print_and_save_log("Output directory: " + output);
print_and_save_log("Input image file suffix: " + suffix);
print_and_save_log("Do manual correction of auto segmentation results?: " + tf(manual_correction));
print_and_save_log("Grid interval (um) and color: " + grid_pars);
print_and_save_log("Min ROI area (um^2) for filtering after manual correction: " + min_roi_area);
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
	open(input + File.separator + file);

	// Load well center ROI if present
	if (isin(sample + "_wellcenter.roi", getFileList(input))) {
		roiManager("Open", input + File.separator + sample + "_wellcenter.roi");
		roiManager("Deselect");
		roiManager("Set Color", "yellow");
		run("From ROI Manager");  // add ROI to overlay
		roiManager("Deselect"); roiManager("Delete");
	}

	// Load auto segmentation ROIs for manual correction
	// Check and ask if continue from last time
	if (isin(out_ROIs_zip, getFileList(output))) {
		if (getBoolean("Detect pre-existing " + out_ROIs_zip + "\n" +
					   "Continue from last time by loading this file?\n" +
					   "Otherwise load " + in_ROIs_zip)) {
			print_and_save_log("Continue from last time ...");
			roiManager("Open", output + File.separator + out_ROIs_zip);
		}
		else roiManager("Open", input + File.separator + in_ROIs_zip);
	}
	else roiManager("Open", input + File.separator + in_ROIs_zip);
	// Print number of ROIs
	print_and_save_log("Initial: " + roiManager("count") + " ROIs");
	// Wait for 1sec to see the image and ROIs if no manual correction and not batch mode hide
	if ((!manual_correction) & (!batch_mode_hide)) {roiManager("show all without labels"); wait(1000);}  
	
	// Manual correction
	if (manual_correction) {
		// Show grid to aid manual correction
		if (grid_pars.length > 0) run("Grid...", "grid=Lines area=" + (pow(parseFloat(grid_pars[0]),2)) + " color=" + (grid_pars[1]) + " center");
		// Manual correction
		roiManager("Show All without labels");
		correct_complete = false;
		while (!correct_complete) {
			waitForUser("Manual correction: \n" +
						"Draw undetected objects: [g] + draw + [t] \n" +
						"Remove incorrect objects: [g] + select + [delete] \n" +
						"Correct detected objects: [g] + draw + [t] + select + add [o] or remove [p] a region \n" +
						"Merge two objects: [g] + select A + [e] + select B + [o] \n" +
						"Split touching objects: [g] or [h] + draw + [t] + select + [u] \n" +
						"Remove small ROIs in the whole image: [q] \n" +
						"Remove small ROIs in an ROI: [g] + draw + [w] \n" +
						"Keyboard shortcut: \n" +
						"- zoom tool [b], scrolling tool [n] \n" +
						"- freehand selection tool [g], freeline selection tool [h] \n" +
						"- ROI manager: show all with labels [j] none [k] all without labels [l] \n" +
						"Click OK when done.");
			correct_complete = getBoolean("Have you completed manual correction?");  // in case of error during manual correction, script proceeds to execute this line, thus user can still go back to continue manual correction and save all previous efforts
		}
		// Print number of ROIs after manual correction
		print_and_save_log("After manual correction: " + roiManager("count") + " ROIs");
		// Filter out manual correction artifacts
		rn = roiManager("count");
		for (r = rn-1; r >= 0; r--) {
		    roiManager("select", r);
			// Mistakenly use freehand line tool to draw an ROI
			if (getValue("Length") > 0) {
				run("To Selection");  // zoom in to selection
				run("Out [-]"); run("Out [-]"); run("Out [-]"); run("Out [-]"); run("Out [-]");
				roiManager("show all with labels");
				setTool("freehand");
				waitForUser("This is a freeline ROI.\n" +
							"Delete it, re-draw by freehand selection tool and add to ROI manager.\n" + 
							"Click OK when done.");
			}
		    // Remove small objects (e.g., due to Split command in "Split ROI" macro)
			else if (getValue("Area") < min_roi_area) roiManager("delete");
		}
		run("Select None"); roiManager("show none");
		// Print number of ROIs after filtering
		print_and_save_log("After removing small ROIs: " + roiManager("count") + " ROIs");
		// Rename ROIs
		rename_rois();
	}
	
	// Remove overlay if there is any
	run("Remove Overlay");
	
	// Save ROIs
	roiManager("Deselect");
	roiManager("Save", output + File.separator + out_ROIs_zip);

	// Remove ROIs on image edges
	rn = roiManager("count");
	for (r = rn-1; r >= 0; r--) {
		roiManager("select", r);
		getSelectionBounds(x, y, w, h);
		if ( (x<=0) || (y<=0) || (x+w>=getWidth()) || (y+h>=getHeight()) ) roiManager("delete");
	}
	run("Select None"); roiManager("show none");
	// Print number of ROIs after manual correction
	print_and_save_log("After excluding ROIs on edges: " + roiManager("count") + " ROIs");
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

/////////////////////////

// Close all windows
function close_everything() {
	// Close all images
	run("Close All"); 
	// Close all non-image windows
	windows = getList("window.titles");
	if (windows.length > 0) {
		for (i = 0; i < windows.length; i++) {
			selectWindow(windows[i]); 
			run("Close");
		}
	}
}

// Print message and save/append it to Log.txt
function print_and_save_log(content) {
	print(content);
	File.append(content, proj + File.separator + "02_Log.txt");	
}

// Strip suffix from the right of a string
function rstrip(string, suffix) {
	return substring(string, 0, lengthOf(string) - lengthOf(suffix));	
}

// Check if a value is in an array
function isin(x, arr) {
	for (i = 0; i < arr.length; i++)
		if (x == arr[i])
			return true; 
	return false; 
}

// Convert true/false to "True"/"False"
function tf(n) { 
	if (n) return "True";
	else return "False";
}

// Rename ROIs in ROI manager to 1-n
function rename_rois() {
	for (r = 0; r < roiManager("count"); r++) {
		roiManager("select", r);
		roiManager("rename", r+1);
	}
	roiManager("show none"); run("Select None");
}