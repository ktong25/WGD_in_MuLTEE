#@ File (label = "Project directory", style = "directory") proj
#@ String (label = "File suffix", value = ".nd2") suffix
#@ Boolean (label = "Output intermediate files", value = false) output_intermediate  // for checking, troubleshooting and demo
#@ Boolean (label = "Use batch mode 'hide'", value = true) batch_mode_hide

// Script info
pipeline = "ploidy_measurement";
version = "v1.2.2";
script = String.join(newArray(pipeline, version), "_");

// Input/Output
input = proj + File.separator + "00_raw";
output = proj + File.separator + "01_output";
if (!File.exists(input)) exit("Input directory 00_raw does not exist.");
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
print_and_save_log("File suffix: " + suffix);
print_and_save_log("Output intermediate files: " + tf(output_intermediate));
print_and_save_log("Use batch mode 'hide': " + tf(batch_mode_hide));
print_and_save_log("");

// Parameters
nucleus_channels = newArray(1,2,3);  // use array, one or more elements
cluster_channel = 4;
min_cluster_area = 500;  //um^2
nucleus_alt_radius = 20;  //px  // auto local threshold radius
min_nucleus_area = 0.5;  //um^2
max_nucleus_area = 10;  //um^2
min_nucleus_circ = 0.7;  // filter out touching nuclei (including most M-phase dividing nuclei) but fine, retain G1/S/G2/early-M round nuclei, eventually use G1 mean
cyto_alt_radius = 20;  //px  // auto local threshold radius

// Initialization
run("Colors...", "foreground=white background=black selection=red");
setBackgroundColor(0, 0, 0);
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
			// Analyze image if this sample has not been completely analyzed (i.e., its final output file does not exist)
			if (!isin(sample + "_hist.png", getFileList(output))) {
				processFile(input, output, file);
			}
		}
	}
}

function processFile(input, output, file) {

	// Print sample info
	print_and_save_log("Processing sample " + sample + " ...");
	
	// Open image
	run("Bio-Formats", "open=" + input + File.separator + file + " autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
	// Check if unit is micron (important for area thresholds etc.)
	getPixelSize(unit, pixelWidth, pixelHeight);
	if (unit != "microns") exit("Unit is not microns.");
	
	// Split channels
	run("Split Channels");
	// If using multiple z slices of nucleus channel, perform max intensity projection
	if (nucleus_channels.length > 1) {
		for (n = 0; n < nucleus_channels.length; n++) {
			selectWindow("C" + nucleus_channels[n] + "-" + file);
			rename("" + (n+1) + "_forstack");
		}
		run("Images to Stack", "  title=forstack");
		run("Z Project...", "projection=[Max Intensity]");
		rename("nucleus_channel");
		if (output_intermediate) save_tiff("nucleus_channel");
		selectWindow("Stack"); close();
	}
	else {
		selectWindow("C" + nucleus_channels[0] + "-" + file);
		rename("nucleus_channel");
		if (output_intermediate) save_tiff("nucleus_channel");
	}
	
	// Cluster segmentation (detect cell region)
	selectWindow("C" + cluster_channel + "-" + file);
	rename("cluster_channel");
	resetMinAndMax(); setOption("ScaleConversions", true); run("8-bit");
	if (output_intermediate) save_tiff("cluster_channel");
	run("Gaussian Blur...", "sigma=0.1 scaled");  //um
	run("Find Edges");
	setAutoThreshold("Li dark");
	setOption("BlackBackground", true);
	run("Convert to Mask");
	run("EDM Binary Operations", "iterations=3 operation=close");
	run("Fill Holes");
	run("Set Measurements...", "area mean shape integrated redirect=None decimal=3");
	run("Analyze Particles...", "size=" + min_cluster_area + "-Infinity show=Masks add in_situ");  // cluster ROIs  // not exclude on edges
	run("Select None"); roiManager("show none"); run("Remove Overlay");
	if (output_intermediate) save_tiff("clusters");
	ncluster = roiManager("count");
	// Report cluster number
	print_and_save_log("Detect " + ncluster + " clusters"); 
	// Rename ROIs as 1-n
	rename_rois();
	// Save cluster ROIs
	roiManager("deselect");
	roiManager("save", output + File.separator + sample + "_cluster_ROIs.zip");
	// Create cell region (all-clusters region)
	roiManager("Select", Array.getSequence(ncluster));
	roiManager("Combine");
	roiManager("Add");  // index of all-clusters ROI is ncluster
	run("Select None"); roiManager("show none");
	// Measure background intensity (method 2) as median intensity of non-cell region
	selectWindow("nucleus_channel");
	roiManager("Select", ncluster);
	run("Make Inverse");
	bgint2 = getValue("Median");  // tolerate few bright small debris
	run("Select None"); roiManager("show none");
	print_and_save_log("Non-cell bgInt: " + bgint2);
	
	// Nucleus segmentation
	// All nuclei (including single or touching nuclei)
	selectWindow("nucleus_channel");
	run("Duplicate...", "title=nuclei_all");
	resetMinAndMax(); setOption("ScaleConversions", true); run("8-bit");
	run("Auto Local Threshold", "method=Phansalkar radius=" + nucleus_alt_radius + " parameter_1=0 parameter_2=0 white");  // px  // auto local threshold is better than auto threshold (on whole image or on one cluster)
	run("Grays");
	run("Analyze Particles...", "size=" + min_nucleus_area + "-Infinity show=Masks in_situ");  // segment all nuclei, touching or not touching // not exclude on edges
	// Clear non-cell region (where there can be foreground pixels like debris or noise due to auto local thresholding)
	roiManager("Select", ncluster);
	run("Clear Outside");
	run("Select None"); roiManager("show none");
	if (output_intermediate) save_tiff("nuclei_all"); 
	// Single nuclei (filtering out touching nuclei)
	run("Duplicate...", "title=nuclei_single");
	run("Analyze Particles...", "size=" + min_nucleus_area + "-" + max_nucleus_area + " circularity=" + min_nucleus_circ + "-1 show=Masks exclude in_situ");  // filter out touching nuclei, exclude on edges
	if (output_intermediate) save_tiff("nuclei_single"); 
	
	// Cytoplasm segmentation (for background subtraction method 1)
	selectWindow("nucleus_channel");
	run("Duplicate...", "title=cyto");
	// Clear nuclei and non-cell region in 16-bit image then convert to 8-bit
	// Why: increase dynamic range of dim cytoplasm pixel values (e.g., from 10-20 to 20-200) to facilitate thresholding
	// Clear nuclei region
	selectWindow("nuclei_all");
	run("Create Selection");
	roiManager("Add");  // index of all-nuclei ROI is ncluster+1
	run("Select None"); roiManager("show none");
	selectWindow("cyto");
	roiManager("select", ncluster+1);
	run("Clear", "slice");
	run("Select None"); roiManager("show none");
	// Clear non-cell region (in case there are bright debris which can affect converting to 8-bit below)
	roiManager("Select", ncluster);
	run("Enlarge...", "enlarge=" + cyto_alt_radius + " pixel");  //px  // clear outside of clusters then auto local threshold will generate some artifact foreground pixels in cluster edges, thus enlarge cluster ROIs and clear outside of enlarged ROIs then auto local threshold and clear outside of original ROIs
	run("Clear Outside");
	run("Select None"); roiManager("show none");
	// Convert to 8-bit with resetting min/max
	run("Enhance Contrast", "saturated=0.35");  // this tolerates few potential super-bright hot spots
	setOption("ScaleConversions", true); run("8-bit");  
	if (output_intermediate) save_tiff("nucleus_channel_cyto");
	// Auto local thresholding
	// Note: there are empty dark space in the crushed clusters, use auto local thresholding to distinguish cytoplasm from them
	run("Auto Local Threshold", "method=Phansalkar radius=" + cyto_alt_radius + " parameter_1=0 parameter_2=0 white");  //px
	run("Grays");
	// Clear non-cell region (also remove cluster-edge foreground pixels which are artifacts)
	roiManager("Select", ncluster);
	run("Clear Outside");
	run("Select None"); roiManager("show none");
	// Clear nuclei region just to make sure
	roiManager("Select", ncluster+1);
	run("Clear", "slice");
	run("Select None"); roiManager("show none");
	if (output_intermediate) save_tiff("cyto");
	
	// For each cluster
	// Measure total fluorescent intensity of each nucleus with background subtraction
	// Measure background intensity as median intensity of (1)all cytoplasm (non-nuclei cell region) in each cluster (2)non-cell region in the whole image
	run("Clear Results");
	c = 0;  // count of nuclei/rows in Results table
	for (r = 0; r < ncluster; r++) {  // for each cluster
		// Measure total fluorescent intensity of each nucleus along with other nucleus metrics and cluster metrics
		selectWindow("nuclei_single");
		roiManager("Select", r);
		run("Set Measurements...", "area mean shape integrated redirect=None decimal=6");
		if (getValue("RawIntDen") > 0) {  // if there are nuclei in the cluster
			cluster_area = getValue("Area");
			cluster_roundness = getValue("Round");
			getSelectionBounds(x, y, w, h);
			if ( (x<=0) || (y<=0) || (x+w>=getWidth()) || (y+h>=getHeight()) ) cluster_onedges = true;
			else cluster_onedges = false;
			run("Set Measurements...", "area mean shape integrated redirect=nucleus_channel decimal=6");  // measure intensity on nucleus_channel image
			run("Analyze Particles...", "display");
			run("Select None"); roiManager("show none");
			updateResults();
			// Measure background intensity (method 1)
			selectWindow("cyto");
			run("Duplicate...", "title=cyto-1");
			roiManager("Select", r);
			run("Clear Outside");
			run("Select None"); roiManager("show none");
			run("Create Selection");
			selectWindow("nucleus_channel");
			run("Restore Selection");
			bgint = getValue("Median");  // background intensity
			run("Select None"); roiManager("show none");
			selectWindow("cyto-1");
			close();
			// For each nucleus, subtract background and add cluster metrics
			for (x = c; x < nResults; x++) {
				area = getResult("Area", x);
				rawarea = round(area / (pixelWidth * pixelHeight));  // thus need more decimal points (6 instead of default 3) in "Set Measurements"
				rawintden = getResult("RawIntDen", x);
				rawintden_bgs = rawintden - rawarea * bgint; 
				rawintden_bgs2 = rawintden - rawarea * bgint2; 
				setResult("cluster_cyto_BackInt", x, bgint);
				setResult("RawIntDen_bgs", x, rawintden_bgs);
				setResult("noncell_BackInt", x, bgint2);
				setResult("RawIntDen_bgs2", x, rawintden_bgs2);
				setResult("cluster_ID", x, r+1);
				setResult("cluster_area", x, cluster_area);
				setResult("cluster_roundness", x, cluster_roundness);
				setResult("cluster_onedges", x, cluster_onedges);
				cluster_nucleinum = nResults - c;
				setResult("cluster_nucleinum", x, cluster_nucleinum);
				updateResults();
			}
			// Report cluster info
			print_and_save_log("Cluster " + (r+1) + ": cyto bgInt " + bgint + " | " + cluster_nucleinum + " nuclei");
			// Update
			c = nResults;
		}
		else {
			print_and_save_log("Cluster " + (r+1) + ": 0 nuclei");
		}
	}
	// Save results
	saveAs("Results", output + File.separator + sample + "_Results.csv");
	
	// Nuclei ROIs
	// Remove pre-existing ROIs
	roiManager("deselect"); roiManager("delete");
	// Get nuclei ROIs
	selectWindow("nuclei_single");
	run("Convert to Mask");
	run("Set Measurements...", "area mean shape integrated redirect=None decimal=3");
	run("Analyze Particles...", "add");
	run("Select None"); roiManager("show none");
	// Report nuclei number
	print_and_save_log("Detect " + roiManager("count") + " nuclei in total"); 
	// Rename ROIs
	rename_rois();
	// Save nucleus ROIs
	roiManager("deselect");
	roiManager("save", output + File.separator + sample + "_nucleus_ROIs.zip");
	
	// Plot nuclei intensity distribution per image
	run("Distribution...", "parameter=RawIntDen_bgs or=50 and=0-0");  // use background subtraction method 1
	saveAs("PNG", output + File.separator + sample + "_hist.png");
	
	// Complete
	// Close all images
	run("Close All"); 
	// Close ROI manager
	roiManager("deselect"); roiManager("delete"); if (!batch_mode_hide) {selectWindow("ROI Manager"); run("Close");}
	// Close Results table
	run("Clear Results"); if (!batch_mode_hide) {selectWindow("Results"); run("Close");}
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
	File.append(content, proj + File.separator + "log.txt");
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

// Convert true/false to "Yes"/"No"
function tf(n) { 
	if (n) return "Yes";
	else return "No";
}

// Rename ROIs in ROI manager to 1-n
function rename_rois() {
	for (r = 0; r < roiManager("count"); r++) {
		roiManager("select", r);
		roiManager("rename", r+1);
	}
	roiManager("show none"); run("Select None");
}

// Save images
function save_tiff(output_suffix) { 
	run("Duplicate...", " ");
	saveAs("TIFF", output + File.separator + sample + "_" + output_suffix + ".tif");
	close();
}