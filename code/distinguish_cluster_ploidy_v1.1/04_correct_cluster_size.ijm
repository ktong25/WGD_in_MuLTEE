#@ File (label = "Project directory", style = "directory") proj

// Script info
pipeline = "distinguish_cluster_ploidy";
version = "v1.1";
step = "04_correct_cluster_size";

// Input/Output
input = proj + File.separator + "01_auto_segmentation";
input2 = proj + File.separator + "02_manual_correction";
input3 = proj + File.separator + "03_measurement";
output = proj + File.separator + "04_measurement_corrected";
if (!File.exists(input)) exit("Input directory 01_auto_segmentation does not exist.");
if (!File.exists(input2)) exit("Input directory 02_manual_correction does not exist.");
if (!File.exists(input3)) exit("Input directory 03_measurement does not exist.");
if (!File.exists(output)) File.makeDirectory(output);

// Parameters
shrink_cluster_radius = 1.8;  //um
suffix = "_downsize.tif";

// Main
setBatchMode("hide");
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
			// Specify input and output file names
			in_ROIs = sample + "_ROIs_auto_manual_eoe.zip";
			out_ROIs = sample + "_ROIs_auto_manual_eoe_corrected.zip";
			in_results = sample + "_Results.csv";
			out_results = sample + "_Results_corrected.csv";
			// Analyze image if this sample has not been completely analyzed (i.e., its final output file does not exist)
			if (!isin(out_results, getFileList(output))) {
				processFile(input, output, file);
			}
		}
	}
}

function processFile(input, output, file) {
	
	// Open input data
	// Image
	open(input + File.separator + file);
	// Check if unit is micron (important for area measurement and shrink_cluster_radius)
	getPixelSize(unit, pixelWidth, pixelHeight);
	if (unit != "microns") exit("Unit is not microns.");
	// ROIs
	roiManager("Open", input2 + File.separator + in_ROIs);
	// Results table
	open(input3 + File.separator + in_results);	
	Table.rename(in_results, "Results");
	
	// Correct cluster size
	// Shrink cluster ROIs to fit the real sizes of clusters, update cluster ROIs and area measurements
	for (r = 0; r < roiManager("count"); r++) {
		roiManager("select", r);
		run("Enlarge...", "enlarge=-" + shrink_cluster_radius);  //um
		roiManager("update");
		area = getValue("Area");
		setResult("Area", r, area);
	}
	run("Select None"); roiManager("show none");
	
	// Save outputs
	// ROIs
	roiManager("Deselect");
	roiManager("Save", output + File.separator + out_ROIs);
	// Results table
	saveAs("Results", output + File.separator + out_results);
	
	// Close
	// Close all images
	run("Close All"); 
	// Close ROI manager
	roiManager("deselect"); roiManager("delete"); 
	// Close Results table
	run("Clear Results"); 

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