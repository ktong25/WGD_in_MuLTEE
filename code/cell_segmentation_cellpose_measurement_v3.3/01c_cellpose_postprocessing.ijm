#@ File (label = "Project directory", style = "directory") proj
#@ String (label = "File suffix", value = "_downsize.tif") suffix
#@ String (label = "Filter cells (um) (e.g., 'Area=2~Inf|Circ.=0.9~1|Skew=-Inf~-1') (if blank then not filter cells)", value = "") cell_filters  // see available measurements in ImageJ below
#@ String (label = "Include only cells within (um) of cluster edges (if blank then not exclude cluster-center cells)", value = "") cluster_edge_width
#@ Boolean (label = "Use batch mode 'hide'", value = true) batch_mode_hide

// Available measurements in ImageJ (see getValue())
// "Area", "Mean", "StdDev", "Mode", "Min", "Max", "X", "Y", "XM", "YM", "Perim.", 
// "BX", "BY", "Width", "Height", "Major", "Minor", "Angle", "Circ.", "Feret", 
// "IntDen", "Median", "Skew", "Kurt", "%Area", "RawIntDen", "Ch", "Slice", "Frame", 
// "FeretX", "FeretY", "FeretAngle", "MinFeret", "AR", "Round", "Solidity", "MinThr", "MaxThr", "Length"

// Script info
pipeline = "cell_segmentation_cellpose_measurement";
version = "v3.3";
step = "01c";
task = "cellpose_postprocessing";
script = String.join(newArray(pipeline, version, step + "_" + task), " ");

// Input/Output
input = proj + File.separator + "01_auto_segmentation";
output = input;
if (!File.exists(input)) exit("Input directory 01_auto_segmentation does not exist.");

// Print and save log
print_and_save_log("--------------------"); 
print_and_save_log("");
print_and_save_log("Running " + script + " ...");
print_and_save_log("ImageJ version: " + getVersion());
print_and_save_log("OS: " + getInfo("os.name"));
print_and_save_log("Project directory: " + proj);
print_and_save_log("Input directory: " + input);
print_and_save_log("Output directory: " + output);
print_and_save_log("File suffix: " + suffix);
print_and_save_log("Filter cells (if blank then not filter cells): " + cell_filters);
print_and_save_log("Include only cells within (um) of cluster edges (if blank then not exclude cluster-center cells): " + cluster_edge_width);
print_and_save_log("Use batch mode 'hide': " + tf(batch_mode_hide));
print_and_save_log("");

// Process input parameter values
// Cell filters: e.g. "Area=2~Inf|Circ.=0.9~1|Skew=-Inf~-1"
cell_filters_split = split(cell_filters, "|");  // e.g. newArray("Area=2~Inf", "Circ.=0.9~1", "Skew=-Inf~-1")
cell_filter_metrics = newArray(0);
cell_filter_mins = newArray(0);
cell_filter_maxs = newArray(0);
for (i = 0; i < cell_filters_split.length; i++) {
	x = split(cell_filters_split[i], "=");  // e.g. newArray("Area", "2~Inf")
	xx = split(x[1], "~");  // e.g. newArray("2", "Inf")
	if ( (xx[0] == "-Inf") && (xx[1] == "Inf") ) exit("Filter cells: min and max cannot simultaneously be -Inf and Inf, respectively");
	cell_filter_metrics = Array.concat(cell_filter_metrics, x[0]);  // e.g. "Area"
	cell_filter_mins = Array.concat(cell_filter_mins, xx[0]);  // e.g. "2"
	cell_filter_maxs = Array.concat(cell_filter_maxs, xx[1]);  // e.g. "Inf"
}
// Cluster edge width
if (cluster_edge_width != "") cluster_edge_width = parseFloat(cluster_edge_width);

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
			// Analyze image if this sample has not been completely analyzed (i.e., its output file does not exist)
			if (!isin(sample + "_masks_auto.png", getFileList(output))) {
				processFile(input, output, file);
			}
		}
	}
}

function processFile(input, output, file) {

	// Print sample info
	print_and_save_log("Processing sample " + sample + " ...");
	
	// Open downsized image
	open(input + File.separator + file);
	rename("Downsize");
	// Check if unit is micron (important for area thresholds, etc.)
	getPixelSize(unit, pixelWidth, pixelHeight);
	if (unit != "microns") exit("Unit is not microns");
	if (pixelWidth != pixelHeight) exit("pixelWidth != pixelHeight");
	
	// Load ROIs
	roiManager("Open", input + File.separator + sample + "_ROIs_cellpose.zip");
	// Report cell number
	ncell_initial = roiManager("count");
	print_and_save_log("Cell number after cellpose segmentation: " + ncell_initial);
	// Rename ROIs
	rename_rois();
	
	// Filter cells using cell-level features
	// For each cell feature filter, count the number of cells that do not pass the filter
	if (cell_filters != "") {
		selectWindow("Downsize");  // filter cells based on .tif (not .png) which keeps metadata like um/px
		cell_filter_counts = newArray(cell_filter_metrics.length); 
		cell_filter_counts = Array.fill(cell_filter_counts, 0);
		rn = roiManager("count");
		for (r = rn-1; r >= 0; r--) {
		    roiManager("select", r);
		    cell_flag = false;
		    for (m = 0; m < cell_filter_metrics.length; m++) {
				metric = cell_filter_metrics[m];
				min = cell_filter_mins[m];
				max = cell_filter_maxs[m];
				// Note: min and max cannot simultaneously be -Inf and Inf, respectively (this is checked at the beginning of this script)
				value = getValue(metric);
				metric_flag = false;
				if ( (min == "-Inf") && (max != "Inf") ) {
					if (value > parseFloat(max)) {
						metric_flag = true; cell_flag = true;
					}
				}
				else if ( (min != "-Inf") && (max == "Inf") ) {
					if (value < parseFloat(min)) {
						metric_flag = true; cell_flag = true;
					}
				}
				else if ( (min != "-Inf") && (max != "Inf") ) {
					if ( (value > parseFloat(max)) || (value < parseFloat(min)) ) {
						metric_flag = true; cell_flag = true;
					}
				}
				if (metric_flag) cell_filter_counts[m]++;
		    }
		    if (cell_flag) roiManager("delete");
		}
		clear_selection();
		// Report cell number
		ncell_filter = roiManager("count");
		print_and_save_log("Cell number after filtering: " + ncell_text(ncell_filter));
		for (m = 0; m < cell_filter_metrics.length; m++) {
			metric = cell_filter_metrics[m]; 
			count = cell_filter_counts[m];
			print_and_save_log("- Cell number not passing " + metric + " filter: " + ncell_text(count));
		}
		// Rename ROIs
		rename_rois();
	}
	
	// Exclude cluster-center cells
	if (cluster_edge_width != "") {
		/// Convert cluster ROIs to binary masks
		selectWindow("Downsize");
		run("Duplicate...", "title=ClusterMask");  // thus um/px info is kept
		run("Select All");
		setBackgroundColor(0, 0, 0); 
		run("Clear", "slice");
		run("Select None");
		rn1 = roiManager("count");
		roiManager("Open", input + File.separator + sample + "_cluster_ROIs.zip");  // append to ROI manager
		rn2 = roiManager("count");
		for (r = rn2-1; r >= rn1; r--) {
		    roiManager("select", r);
		    setForegroundColor(255, 255, 255);
			run("Fill", "slice");
			roiManager("delete");
		}
		clear_selection();
		/// Erode cluster masks and subtract to get cluster edges
		/// Note: somehow BioVoxxel's EDM Binary Operations cannot perform more than 254 cycles of erosion at one time, thus do it several times, each time erode 200 cycles
		/// Note: it is possible that there is nothing left after erosion but this is okay
		run("Duplicate...", "title=ClusterMaskErode");
		getPixelSize(unit, pixelWidth, pixelHeight);
		cluster_edge_erode_cycles = round(cluster_edge_width / pixelWidth); 
		run("Options...", "iterations=1 count=1 black pad do=Nothing");  // pad edges when eroding thus do not erode the side of cluster that touches image edge
		remaining_cycles = cluster_edge_erode_cycles;
		while (remaining_cycles > 0) {
			cur_cycles = minOf(remaining_cycles, 200);
			run("EDM Binary Operations", "iterations=" + cur_cycles + " operation=erode");
			remaining_cycles = remaining_cycles - 200;
		}
		imageCalculator("Subtract", "ClusterMask", "ClusterMaskErode");  // directly modify original image
		selectWindow("ClusterMaskErode"); close();
		selectWindow("ClusterMask"); rename("ClusterEdgeMask");
		/// Include only cluster-edge cells (cell ROIs overlapping with (but not necessarily entirely within) the cluster-edge masks)
		rn = roiManager("count");
		for (r = rn-1; r >= 0; r--) {
			roiManager("select", r);
			if (getValue("%Area") == 0) roiManager("delete");
		}
		clear_selection();
		/// Close temporary images
		selectWindow("ClusterEdgeMask"); close();
		// Report cell number
		ncell_edge = roiManager("count");
		print_and_save_log("Cell number after excluding cluster-center cells: " + ncell_text(ncell_edge));
		// Rename ROIs
		rename_rois();
	}

	// Save cell ROIs
	roiManager("deselect");
	roiManager("Save", output + File.separator + sample + "_ROIs_auto.zip");
	
	// Convert cell ROIs to label image and save
	selectWindow("Downsize");
	roiManager("deselect");
	run("ROI Map");  // create a new image named "Roi Map"
	run("Grays");
	saveAs("PNG", output + File.separator + sample + "_masks_auto.png");
	close();
	
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

// Print summary text for cell number
function ncell_text(ncell) {
	return "" + ncell + " (" + percentage(ncell, ncell_initial, 1) + "%)";  // percentage relative to initial cell number
}
