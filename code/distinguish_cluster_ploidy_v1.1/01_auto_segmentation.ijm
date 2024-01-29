#@ File (label = "Project directory", style = "directory") proj
#@ String (label = "File suffix", value = ".nd2") suffix
#@ Boolean (label = "Remove numerical prefix of input file names when output (e.g., 01_t0 -> t0)", value = true) remove_num_prefix
#@ Float (label = "Min cluster area (um^2)", value = 50) min_cluster_area
#@ Boolean (label = "Use batch mode 'hide'", value = true) batch_mode_hide

// Script info
pipeline = "distinguish_cluster_ploidy";
version = "v1.1";
step = "01_auto_segmentation";
script = String.join(newArray(pipeline, version, step), "_");

// Input/Output
input = proj + File.separator + "00_raw";
output = proj + File.separator + "01_auto_segmentation";
if (!File.exists(input)) exit("Input directory 00_raw does not exist.");
if (!File.exists(output)) File.makeDirectory(output);

// Print and save log
print_and_save_log("--------------------"); 
print_and_save_log("");
print_and_save_log("Running " + script + " ...");
print_and_save_log("ImageJ version: " + getVersion());
print_and_save_log("Project directory: " + proj);
print_and_save_log("Input directory: " + input);
print_and_save_log("Output directory: " + output);
print_and_save_log("File suffix: " + suffix);
print_and_save_log("Remove numerical prefix of input file names when output: " + tf(remove_num_prefix));
print_and_save_log("Min cluster area (um^2): " + min_cluster_area);
print_and_save_log("Use batch mode 'hide': " + tf(batch_mode_hide));
print_and_save_log("");

// Initialization
run("Colors...", "foreground=white background=black selection=red");
run("Options...", "iterations=1 count=1 black pad");  // binary operations
if (batch_mode_hide) setBatchMode("hide");  // set batch mode

// Main
close_everything();
p = -1; 
processFolder(input, p);
print_and_save_log("--------------------"); 
close_everything();

function processFolder(input, p) {
	p++;
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		file = list[i];
		if (File.isDirectory(input + File.separator + file))
			processFolder(input + File.separator + file, p);
		if (endsWith(file, suffix)) {
			// Process sample name
			/// Create prefix
			os = getInfo("os.name");
			if (startsWith(os, "Windows")) input_split = split(rstrip(input, "/"), "\\");
			else if (startsWith(os, "Mac OS")) input_split = split(input, "/");
			prefix = String.join(Array.slice(input_split, input_split.length-p, input_split.length), "_");
			if (p != 0) prefix += "_";
			/// Extract file name
			filename = rstrip(file, suffix);  // strip file extension
			if (remove_num_prefix) filename = substring(filename, indexOf(filename, "_")+1);  // strip numerical prefix, e.g. 01_t0 -> t0; if not, t0 -> t0, 01_t0 -> 01_t0
			/// Combine
			sample = prefix + filename;
			// Analyze image if this sample has not been completely analyzed (i.e., its output file does not exist)
			if (!isin(sample + "_ROIs_auto.zip", getFileList(output))) {
				processFile(input, output, file);
			}
		}
	}
}

function processFile(input, output, file) {

	// Print sample info
	print_and_save_log("Processing sample " + sample + " ...");

	// Open raw image
	run("Bio-Formats", "open=" + input + File.separator + file + " autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
	rename(sample);
	// Check if unit is micron (important for area thresholds, etc.)
	getPixelSize(unit, pixelWidth, pixelHeight);
	if (unit != "microns") exit("Unit is not microns.");
	
	// Stack to images
	getDimensions(width, height, nchannel, nslice, nframe);
	if (nslice == 1) rename(sample + "-0001");
	else run("Stack to Images");
	// Select z slice for cluster segmentation
	middle_slice = round(nslice/2); // 1->1, 2->1, 3->2, 4->2
	selectWindow(sample + "-000" + middle_slice);
	
	// Cluster segmentation
	// Convert to 8-bit
	resetMinAndMax(); setOption("ScaleConversions", true); run("8-bit");
	// Save downsized image for manual correction
	saveAs("TIFF", output + File.separator + sample + "_downsize.tif");  // window name changes to XXX.tif
	// Segmentation
	run("Gaussian Blur...", "sigma=1");  //px
	run("Find Edges");
	run("Auto Local Threshold", "method=Bernsen radius=5 parameter_1=0 parameter_2=0 white");  //px
	run("EDM Binary Operations", "iterations=4 operation=close");
	run("Fill Holes");
	run("Analyze Particles...", "size=" + min_cluster_area + "-Infinity add");  // cluster ROIs are slightly larger than clusters but fine
	roiManager("show none"); run("Select None");

	// Print number of ROIs after filtering
	print_and_save_log("Auto segmentation: " + roiManager("count") + " ROIs");
	// Rename ROIs
	rename_rois();
	// Save cluster ROIs
	roiManager("Deselect");
	roiManager("Save", output + File.separator + sample + "_ROIs_auto.zip");
		
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
	File.append(content, proj + File.separator + "01_Log.txt");	
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