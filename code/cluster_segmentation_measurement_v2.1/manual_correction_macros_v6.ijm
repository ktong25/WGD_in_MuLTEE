// Below is a collection of macros for aiding manual correction of segmentation
// Install: copy and paste all the following code into "StartupMacros.fiji.ijm" file in "Fiji.app/macros" folder (show by File -> Show Folder -> Macros), restart ImageJ
// Note: these may override ImageJ's default keyboard shortcuts, e.g. [n] is used to create a new image in ImageJ by default

// Clear all selections
function clear_selection() { 
	roiManager("deselect"); run("Select None"); roiManager("show none");
}

function an_ROI_selected() { 
	if (roiManager("index") != -1) return true;
	else return false;
}

function a_non_ROI_selection() { 
	if ( (selectionType() != -1) & (roiManager("index") == -1) ) return true;
	else return false;
}

macro "Delete ROI [x]" {  // select an ROI and hit keyboard shortcut
	if (an_ROI_selected()) {  // if one ROI is selected, delete it; otherwise do nothing (because if no ROI is selected, delete will be deleting all ROIs)
		roiManager("delete");
		roiManager("show all with labels");
	}
}

// Macro "Split ROI" with keyboard shortcut [u]
// Last modified: Kai Tong, 2022/3/4
// Tested ImageJ version: 2.3.0/1.53f
// Function: split an ROI into two or more ROIs using a manually-drawn split line or shape
// Install: copy and paste the following code into "StartupMacros.fiji.ijm" file in "Fiji.app/macros" folder (show by File -> Show Folder -> Macros), restart ImageJ
// Usage: 
// - open an image and its ROIs-containing ROI manager
// - find an ROI you want to split
// - use freehand line tool or freehand selection tool to draw a split line or shape (both can go beyond the focal ROI), add it to ROI manager (shortcut [t])
// - use freehand line/selection tool to select the ROI to be split (when "Show All with labels" in ROI manager is on, click the number at the center of the ROI in the image), or select it in ROI manager
// - press shortcut key to run the macro
// Algorithm: perform AND, XOR, Split to get post-split ROIs, delete original and intermediate ROIs
// Notes: 
// - customize keyboard shortcut below where indicated
// - using split line can split one ROI into multiple ROIs, but post-split ROIs may still be "connected" and may be recognized as one object by Analyze Particles

macro "Split ROI [u]" {
	if (an_ROI_selected()) {
		a = roiManager("index");  // index of ROI (selected) to be split
		b = roiManager("count") - 1;  // index of freehand line/selection ROI that is used to split (previously added to the end of ROI manager)
		roiManager("Select", newArray(a,b));
		roiManager("AND");
		roiManager("Add");
		c = b + 1; 
		roiManager("Select", newArray(a,c));
		roiManager("XOR");
		roiManager("Add");
		d = c + 1; 
		roiManager("Select", d);
		roiManager("Split");
		roiManager("Select", newArray(a,b,c,d));
		roiManager("Delete");
		roiManager("Show All with labels");
	}
}

// Macro "Add to ROI" with keyboard shortcut [o]
// Last modified: Kai Tong, 2022/3/4
// Tested ImageJ version: 2.3.0/1.53f
// Function: add a region to an ROI using a manually-drawn shape
// Install: copy and paste the following code into "StartupMacros.fiji.ijm" file in "Fiji.app/macros" folder (show by File -> Show Folder -> Macros), restart ImageJ
// Usage: 
// - open an image and its ROIs-containing ROI manager
// - find an ROI you want to add to
// - use freehand selection tool to draw a shape (it must partially overlap with the focal ROI, and the region that is not overlapping is used to add to the focal ROI), add it to ROI manager (shortcut [t])
// - use freehand selection tool to select the ROI to be added to (when "Show All with labels" in ROI manager is on, click the number at the center of the ROI in the image), or select it in ROI manager
// - press shortcut key to run the macro
// Algorithm: perform OR to get the result ROI, delete original and intermediate ROIs
// Notes: 
// - customize keyboard shortcut below where indicated

macro "Add to ROI [o]" {
	if (an_ROI_selected()) {
		a = roiManager("index");  // index of ROI (selected) to be added to
		b = roiManager("count") - 1;  // index of freehand selection ROI that is used to add (previously added to the end of ROI manager)
		roiManager("Select", newArray(a,b));
		roiManager("Combine");
		roiManager("Add");
		roiManager("Select", newArray(a,b));
		roiManager("Delete");
		roiManager("Show All with labels");
	}
}

// Macro "Remove from ROI" with keyboard shortcut [p]
// Last modified: Kai Tong, 2022/3/4
// Tested ImageJ version: 2.3.0/1.53f
// Function: remove a region from an ROI using a manually-drawn shape
// Install: copy and paste the following code into "StartupMacros.fiji.ijm" file in "Fiji.app/macros" folder (show by File -> Show Folder -> Macros), restart ImageJ
// Usage: 
// - open an image and its ROIs-containing ROI manager
// - find an ROI you want to remove from
// - use freehand selection tool to draw a shape (it must partially overlap with the focal ROI, and the overlapping region is used to remove from the focal ROI), add it to ROI manager (shortcut [t])
// - use freehand selection tool to select the ROI to be removed from (when "Show All with labels" in ROI manager is on, click the number at the center of the ROI in the image), or select it in ROI manager
// - press shortcut key to run the macro
// Algorithm: perform AND, XOR to get the result ROI, delete original and intermediate ROIs
// Notes: 
// - customize keyboard shortcut below where indicated

macro "Remove from ROI [p]" {
	if (an_ROI_selected()) {
		a = roiManager("index");  // index of ROI (selected) to be removed from
		b = roiManager("count") - 1;  // index of freehand selection ROI that is used to remove (previously added to the end of ROI manager)
		roiManager("Select", newArray(a,b));
		roiManager("AND");
		roiManager("Add");
		c = b + 1; 
		roiManager("Select", newArray(a,c));
		roiManager("XOR");
		roiManager("Add");
		roiManager("Select", newArray(a,b,c));
		roiManager("Delete");
		roiManager("Show All with labels");
	}
}

macro "Move ROI to the end of ROI manager [e]" {   // select an ROI and hit keyboard shortcut
	// Tip: this is helpful for merging two ROIs, by [g] + select + [e] + select + [o]
	if (an_ROI_selected()) {
		a = roiManager("index");
		roiManager("Select", a);
		roiManager("Add");
		roiManager("Select", a);
		roiManager("Delete");
		roiManager("Show All with labels");
	}
}

macro "Merge two ROIs by filling 1px gap [y]" {   // [g] + select + [e] + select + [y]
	// Tip: fill 1px-line gap, typically used to merge two ROIs that come from an incorrectly-split object
	// Note: if the gap is more than 1px, then this macro basically does nothing (besides the slight change due to close operation)
	if (an_ROI_selected()) {
		a = roiManager("index");  // index of ROI selected
		b = roiManager("count") - 1;  // index of another ROI that is to be merged (previously moved to the end of ROI manager)
		roiManager("select", newArray(a,b));
		roiManager("combine");
		// Duplicate two-ROI region
		image = getTitle();
		run("Duplicate...", " title=1");
		run("Duplicate...", " title=2");
		// Raw image
		selectWindow("1");
		clear_selection();
		// Mask image
		selectWindow("2");
		setForegroundColor(255, 255, 255);
		run("Fill", "slice");
		setBackgroundColor(0, 0, 0);
		run("Clear Outside");
		clear_selection();
		getDimensions(width, height, channels, slices, frames);
		extension = 3;  // px 
		run("Canvas Size...", "width=" + (width+extension*2) + " height=" + (height+extension*2) + " position=Center zero");  // this is to avoid image edge artifacts caused by close operation
		run("EDM Binary Operations", "iterations=1 operation=close");  // fill 1px-line gap between two ROIs
		makeRectangle(extension, extension, width, height); run("Crop");
		// Generate ROI based on original image
		selectWindow("2"); run("Copy"); selectWindow(image); run("To Bounding Box"); run("Paste");
		run("Analyze Particles...", "size=0-Infinity add");
		selectWindow("1"); run("Copy"); selectWindow(image); run("Paste");
		clear_selection();
		// Close
		selectWindow("1"); close();
		selectWindow("2"); close();
		// Delete previous ROIs
		roiManager("select", newArray(a,b));
		roiManager("delete");
		// Show ROIs
		selectWindow(image);
		roiManager("show all with labels");
	}
}

macro "Delete small ROIs in the image [q]" {
	// Tip: remove manual correction artifacts (e.g., delete few-pixel objects due to Split command in "Split ROI" macro)
	min_roi_area = 40;  // um^2  // adjust if needed
	rn = roiManager("count");
	for (r = rn-1; r >= 0; r--) {
	    roiManager("select", r);
		if (getValue("Area") < min_roi_area) roiManager("delete");
	}
	clear_selection(); 
	roiManager("Show All with labels");
}

macro "Delete all ROIs in a selection [w]" {  // create a selection (e.g., rectangle, freehand) and hit keyboard shortcut
	// Note: The ROIs to be deleted should be entirely inside the selection
	// Algorithm: 
	// - For all ROIs within a selection, their centroids must be inside the bounding box of the selection
	// - This can be used as an initial, faster screening, before the subsequent, slower final determination using ROI combine
	if (a_non_ROI_selection()) {
		roiManager("add");
		draw_roi_area = getValue("Area");
		getSelectionBounds(X, Y, W, H); X1 = X; Y1= Y; X2 = X+W; Y2 = Y+H; 
		rn = roiManager("count") - 1;  // exclude the user-drawn ROI
		for (r = rn-1; r >= 0; r--) {
		    roiManager("select", r);
		    getSelectionBounds(x, y, w, h); x0 = round(x+w/2); y0 = round(y+h/2);
			if ( (X1<=x0) & (x0<=X2) & (Y1<=y0) & (y0<=Y2) ) {
				roiManager("select", newArray(r, roiManager("count") - 1));
				roiManager("combine");
				if (getValue("Area") == draw_roi_area) {clear_selection(); roiManager("select", r); roiManager("delete");}
			}
		}
		clear_selection();
		roiManager("select", roiManager("count") - 1); roiManager("delete");
		clear_selection();
		roiManager("Show All with labels");
	}
}

// Below set custom keyboard shortcuts for several ImageJ commands to aid manual correction. 

macro "Activate zoom tool [b]" {
	setTool("zoom");
}

macro "Activate hand tool [n]" {  // also called scrolling tool
	setTool("hand");
}

macro "Activate rectangle tool [f]" {
	setTool("rectangle");
}

macro "Activate freehand selection tool [g]" {
	setTool("freehand");
}

macro "Activate freehand line tool [h]" {
	setTool("freeline");
}

macro "Show All with labels [j]" {
	roiManager("Show All with labels");
}

macro "Show None [k]" {
	roiManager("Show None");
}

macro "Show All without labels [l]" {
	roiManager("Show All without labels");
}

// Close all windows
macro "Close all windows [L]" {
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

// Below are manual correction macros specifically designed for aiding cluster segmentation

macro "Auto segment cluster(s) in a selection [r]" {  // create a selection (e.g., rectangle, freehand) (leave some background area and shading should not be too variable) and hit keyboard shortcut
	if (a_non_ROI_selection()) {
		// Parameters
		edm_close_cycles = 4;
		min_cluster_area = 40;  //um^2
		// Main
		image = getTitle();
		run("Duplicate...", "title=1");  // Raw
		run("Duplicate...", "title=2");  // Segmentation
		// Segmentation using auto thresholding (thus some background in the user-drawn ROI is necessary and shading should not be too variable)
		selectWindow("2"); 
		setAutoThreshold("Default");
		setOption("BlackBackground", true);
		run("Convert to Mask");
		run("EDM Binary Operations", "iterations=" + edm_close_cycles + " operation=close");
		run("Fill Holes");
		// Generate ROI based on original image
		selectWindow("2"); run("Copy"); selectWindow(image); run("Paste");
		run("Analyze Particles...", "size=" + min_cluster_area + "-Infinity add");
		selectWindow("1"); run("Copy"); selectWindow(image); run("Paste");
		// Close
		selectWindow("1"); close();
		selectWindow("2"); close();
		// Show ROIs
		selectWindow(image);
		roiManager("show all with labels");
	}
}

macro "Auto split touching clusters [i]" {  // select an ROI and hit keyboard shortcut
	if (an_ROI_selected()) {
		// Parameters
		gaussian_filter_radius = 7;  //um
		find_maxima_prominence = 5;
		min_cluster_area = 40;  //um^2
		// Main
		r = roiManager("index");
		image = getTitle();
		run("Duplicate...", " title=1");
		run("Duplicate...", " title=2");
		run("Duplicate...", " title=3");
		// Raw image
		selectWindow("1");
		clear_selection();
		// Marker image
		selectWindow("2");
		clear_selection();
		run("Gaussian Blur...", "sigma=" + gaussian_filter_radius + " scaled");  //um
		run("Find Maxima...", "prominence=" + find_maxima_prominence + " strict exclude light output=[Single Points]");  // generate "2 Maxima" image window
		// Mask image
		selectWindow("3");
		setForegroundColor(255, 255, 255);
		run("Fill", "slice");
		setBackgroundColor(0, 0, 0);
		run("Clear Outside");
		clear_selection();
		// Marker-controlled watershed
		run("Marker-controlled Watershed", "input=1 marker=[2 Maxima] mask=3 compactness=0 binary calculate use"); // generate "1-watershed" image window (32-bit, where background is 0, objects are labeled 1-n)
		setOption("ScaleConversions", false); run("8-bit"); setOption("ScaleConversions", true);
		setThreshold(1, 255, "raw"); setOption("BlackBackground", true); run("Convert to Mask");
		// Generate ROI based on original image
		selectWindow("1-watershed"); run("Copy"); selectWindow(image); run("Paste");
		run("Analyze Particles...", "size=" + min_cluster_area + "-Infinity add");
		selectWindow("1"); run("Copy"); selectWindow(image); run("Paste");
		clear_selection();
		// Close
		selectWindow("1"); close();
		selectWindow("2"); close();
		selectWindow("3"); close();
		selectWindow("2 Maxima"); close();
		selectWindow("1-watershed"); close();
		// Delete previous ROI
		roiManager("select", r); roiManager("delete");
		// Show ROIs
		selectWindow(image);
		roiManager("show all with labels");
	}
}

// Below are manual correction macros specifically designed for correcting cluster segmentation with many debris detected

// Cluster and debris are assigned and distinguished by ROI color
// ROI color: cluster "red" (same as the default selection color in "ImageJ -> Edit -> Options -> Colors..."), debris "white"
// All white ROIs are: initially auto assigned as white, or later manually assigned as white
// Cluster ROI colors can be red (initially auto or later manually assigned as red) or none (by default, including any newly added ROIs)

macro "Classify all ROIs in the image as cluster/debris by area [Q]" {  // hit keyboard shortcut
	// Algorithm: ROIs below certain area threshold are preliminarily assigned as debris
	// Parameters
	large_roi_area_threshold = 1500;  //um^2
	// Main
	for (r = 0; r < roiManager("count"); r++) {
		roiManager("select", r);
		if (getValue("Area") < large_roi_area_threshold) roiManager("Set Color", "white");
	}
	clear_selection();
	roiManager("show all without labels");  // thus see ROI color more easily
}

// Below two macros: select an ROI or create a selection (e.g., rectangle, freehand) and hit keyboard shortcut
// Select an ROI: assign that ROI as a certain class
// Create a selection: assign all ROIs entirely within that selection as a certain class
// Assign ROI(s) as a certain class: set ROI color

macro "Assign ROI(s) as cluster [1]" {
	assign_roi_class("red");
}

macro "Assign ROI(s) as debris [2]" { 
	assign_roi_class("white");
}

function assign_roi_class(roi_color) { 
	if (an_ROI_selected()) {
		roiManager("Set Color", roi_color);
	}
	else if (a_non_ROI_selection()) {
		roiManager("add");
		draw_roi_area = getValue("Area");
		getSelectionBounds(X, Y, W, H); X1 = X; Y1= Y; X2 = X+W; Y2 = Y+H; 
		rn = roiManager("count") - 1;  // exclude the user-drawn ROI
		for (r = rn-1; r >= 0; r--) {
		    roiManager("select", r);
		    getSelectionBounds(x, y, w, h); x0 = round(x+w/2); y0 = round(y+h/2);
			if ( (X1<=x0) & (x0<=X2) & (Y1<=y0) & (y0<=Y2) ) {
				roiManager("select", newArray(r, roiManager("count") - 1));
				roiManager("combine");
				if (getValue("Area") == draw_roi_area) {clear_selection(); roiManager("select", r); roiManager("Set Color", roi_color); clear_selection();}
			}
		}
		clear_selection();
		roiManager("select", roiManager("count") - 1); roiManager("delete");  // delete the user-drawn ROI
	}
	clear_selection();
	roiManager("show all without labels");  // thus see ROI color more easily
}

macro "Delete all debris ROIs in the image [Z]" {  // hit keyboard shortcut
	rn = roiManager("count");
	for (r = rn-1; r >= 0; r--) {
		roiManager("select", r);
		if (Roi.getStrokeColor == "white") roiManager("delete");
	}
	clear_selection();
	roiManager("show all without labels");  // thus see ROI color more easily
}
