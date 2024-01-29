#@ File (label = "Project directory", style = "directory") proj
#@ String (label = "File suffix", value = ".nd2") suffix
#@ Boolean (label = "Remove numerical prefix of input file names when output (e.g., 01_t0 -> t0)", value = true) remove_num_prefix
#@ Boolean (label = "Use batch mode 'hide'", value = true) batch_mode_hide

// Script info
pipeline = "distinguish_cluster_ploidy";
version = "v1.1";
step = "03_measurement";
script = String.join(newArray(pipeline, version, step), "_");

// Input/Output
input = proj + File.separator + "00_raw";
input2 = proj + File.separator + "02_manual_correction";
output = proj + File.separator + "03_measurement";
if (!File.exists(input)) exit("Input directory 00_raw does not exist.");
if (!File.exists(input2)) exit("Input directory 02_manual_correction does not exist.");
if (!File.exists(output)) File.makeDirectory(output);

// Print and save log
print_and_save_log("--------------------"); 
print_and_save_log("");
print_and_save_log("Running " + script + " ...");
print_and_save_log("ImageJ version: " + getVersion());
print_and_save_log("Project directory: " + proj);
print_and_save_log("Input directory for raw images: " + input);
print_and_save_log("Input directory for ROIs: " + input2);
print_and_save_log("Output directory: " + output);
print_and_save_log("File suffix: " + suffix);
print_and_save_log("Remove numerical prefix of input file names when output: " + tf(remove_num_prefix));
print_and_save_log("Use batch mode 'hide': " + tf(batch_mode_hide));
print_and_save_log("Area unit: um^2");
print_and_save_log("");

// Parameters
min_cluster_area = 50;  //um^2
min_cell_area = 2;  //um^2
max_cell_area = 40;  //um^2
min_cell_circ = 0.7; 
cluster_edge_radius = 10;  //um

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
			// Specify input and output file names
			in_ROIs_eoe_zip = sample + "_ROIs_auto_manual_eoe.zip";
			out_results = sample + "_Results.csv";
			out_plot = sample + "_Plot.png";
			// Analyze image if this sample has not been completely analyzed (i.e., its final output file does not exist)
			if (!isin(out_plot, getFileList(output))) {
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
	// Check if unit is micron (important for area measurement, etc.)
	getPixelSize(unit, pixelWidth, pixelHeight);
	if (unit != "microns") exit("Unit is not microns.");
	
	// Load cluster ROIs
	roiManager("Open", input2 + File.separator + in_ROIs_eoe_zip);
	// Print number of ROIs
	print_and_save_log("" + roiManager("count") + " ROIs");
	
	// Cluster-edge cell segmentation
	// Convert to 8-bit
	selectWindow(sample);
	resetMinAndMax(); setOption("ScaleConversions", true); run("8-bit");
	// Stack to images
	getDimensions(width, height, nchannel, nslice, nframe);
	if (nslice == 1) rename(sample + "-0001");
	else run("Stack to Images");
	// Process each z slice
	for (z = 1; z <= nslice; z++) {
		image = sample + "-000" + z;
		selectWindow(image);
		// Cell segmentation
		// Note: loose filtering + fill holes + strict filtering allows detecting cells with large vacuoles (which is initially objects with holes thus lower circ and smaller area)
		// Note: running watershed or close operation is risky here, thus not aim to split doublets or fix broken cells
		run("Auto Local Threshold", "method=Bernsen radius=20 parameter_1=0 parameter_2=0 white");  //px
		run("Analyze Particles...", "size=0-" + max_cell_area + " circularity=0-1.00 show=Masks exclude in_situ"); 
		run("Fill Holes");
		run("Analyze Particles...", "size=" + min_cell_area + "-" + max_cell_area + " circularity=" + min_cell_circ + "-1.00 show=Masks exclude in_situ");
		// Clear outside of clusters
		roiManager("Select", Array.getSequence(roiManager("count")));
		roiManager("Combine");  // create all-clusters selection
		setBackgroundColor(0, 0, 0);
		run("Clear Outside");
		run("Select None");
		// Filter cluster-edge cells
		/// Create cluster edge selection of all clusters
		selectWindow(image);
		run("Duplicate...", " "); 
		run("Restore Selection");  // restore all-clusters selection
		setForegroundColor(255, 255, 255);
		run("Fill", "slice");
		run("Enlarge...", "enlarge=-" + cluster_edge_radius);  //um
		run("Clear", "slice");
		run("Select None");
		run("Create Selection");  // create all-clusters edge selection
		run("Select None");
		close();
		/// Create image containing cluster-center cells
		selectWindow(image);
		run("Duplicate...", " ");
		run("Restore Selection");  // restore all-clusters edge selection
		run("Fill", "slice");  // this makes cluster edges and all cells touching/within cluster edges white
		run("Select None");
		run("Analyze Particles...", "size=" + min_cell_area + "-" + max_cell_area + " circularity=" + min_cell_circ + "-1.00 show=Masks exclude in_situ");  // filter out cluster edges which should be much larger than cells and have irregular shapes
		/// Subtract image to filter out cluster center cells
		imageCalculator("Subtract", image, image + "-1");
		selectWindow(image + "-1");
		close();
		// Save
		saveAs("Tiff", output + File.separator + sample + "_z" + z + "_cells.tif");  // change image window name to XXX.tif
	}
		
	// Measure cell size in each cluster
	/// Measure
	run("Set Measurements...", "area redirect=None decimal=3");
	ncells = newArray(0);
	cell_mean_areas = newArray(0);
	cell_top5_mean_areas = newArray(0);
	rn = roiManager("count");  // cluster number
	run("Clear Results");
	for (r = 0; r < rn; r++) {
		for (z = 1; z <= nslice; z++) {
			image = sample + "_z" + z + "_cells.tif";
			selectWindow(image);
			roiManager("select", r);
			run("Analyze Particles...", "display");
		}
		ncell = nResults;
		if (ncell > 0) {
			cell_areas = Table.getColumn("Area", "Results");
			Array.getStatistics(cell_areas, min, max, cell_mean_area, stdDev);
			cell_top5_areas = Array.trim(Array.reverse(Array.sort(cell_areas)), 5);  // if input array has fewer than 5 elements then output array has fewer than 5 elements
			Array.getStatistics(cell_top5_areas, min, max, cell_top5_mean_area, stdDev);
			ncells = Array.concat(ncells, ncell);
			cell_mean_areas = Array.concat(cell_mean_areas, cell_mean_area);
			cell_top5_mean_areas = Array.concat(cell_top5_mean_areas, cell_top5_mean_area);
		}
		else {
			ncells = Array.concat(ncells, ncell);
			cell_mean_areas = Array.concat(cell_mean_areas, 0);
			cell_top5_mean_areas = Array.concat(cell_top5_mean_areas, 0);
		}
		run("Clear Results");
	}
	run("Select None"); roiManager("show none");
	/// Fill results table with measurements
	roiManager("deselect");
	run("Set Measurements...", "area redirect=None decimal=3");
	roiManager("measure");
	for (r = 0; r < rn; r++) {
		setResult("cluster_id", r, r+1);
		setResult("ncell", r, ncells[r]);
		setResult("cell_mean_area", r, cell_mean_areas[r]);
		setResult("cell_top5_mean_area", r, cell_top5_mean_areas[r]);
	}
	/// Save results
	if (!batch_mode_hide) updateResults();
	saveAs("Results", output + File.separator + out_results);
	
	// Make plot
	Plot.create("Plot of Results", "Area", "cell_top5_mean_area");
	Plot.add("Circle", Table.getColumn("Area", "Results"), Table.getColumn("cell_top5_mean_area", "Results"));
	Plot.setStyle(0, "blue,#a0a0ff,1.0,Circle");
	Plot.show();
	saveAs("PNG", output + File.separator + out_plot);
	
	// Close
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
	File.append(content, proj + File.separator + "03_Log.txt");	
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

// Get index of a value in an array
function index(x, arr) { 
	for (i = 0; i < arr.length; i++) {
		if (x == arr[i]) {
			return i;
		}
	}
}