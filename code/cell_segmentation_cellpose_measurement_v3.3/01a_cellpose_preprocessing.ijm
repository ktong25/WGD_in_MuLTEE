#@ File (label = "Project directory", style = "directory") proj
#@ String (label = "File suffix", value = ".nd2") suffix
#@ Boolean (label = "Remove numerical prefix of input file names when output (e.g., 01_t0 -> t0)", value = false) remove_num_prefix
#@ String (label = "Position of channel for cell segmentation", value = "1") seg_channel_pos
#@ String (label = "Position of BF channel for cluster segmentation (if blank then skip cluster segmentation)", value = "") cluster_channel_pos  // if "" then do not do this
#@ Boolean (label = "Clear non-cluster region in the channel for cell segmentation", value = false) clear_non_cluster
#@ Boolean (label = "Crop images", value = false) crop_images
#@ Boolean (label = "Use batch mode 'hide' (not if crop images)", value = true) batch_mode_hide

// Check input parameter values
if (crop_images && batch_mode_hide) exit("Cannot use batch mode 'hide' if crop images.");

// Script info
pipeline = "cell_segmentation_cellpose_measurement";
version = "v3.3";
step = "01a";
task = "cellpose_preprocessing";
script = String.join(newArray(pipeline, version, step + "_" + task), " ");

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
print_and_save_log("OS: " + getInfo("os.name"));
print_and_save_log("Project directory: " + proj);
print_and_save_log("Input directory: " + input);
print_and_save_log("Output directory: " + output);
print_and_save_log("File suffix: " + suffix);
print_and_save_log("Remove numerical prefix of input file names when output: " + tf(remove_num_prefix));
print_and_save_log("Position of channel for cell segmentation: " + seg_channel_pos);
print_and_save_log("Position of BF channel for cluster segmentation (if blank then skip cluster segmentation): " + cluster_channel_pos);
print_and_save_log("Clear non-cluster region in the channel for cell segmentation: " + tf(clear_non_cluster));
print_and_save_log("Crop images: " + tf(crop_images));
print_and_save_log("Use batch mode 'hide' (not if crop images): " + tf(batch_mode_hide));
print_and_save_log("");

// Parameters
min_cluster_area = 500;  //um^2
max_cluster_edgehole_area = 150;  //um^2

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
		if (startsWith(getInfo("os.name"), "Windows")) if (endsWith(file, "/")) file = rstrip(file, "/");  // if Windows directory
		if (File.isDirectory(input + File.separator + file))
			processFolder(input + File.separator + file, p);
		if (endsWith(file, suffix)) {
			// Process sample name
			/// Create prefix
			input_split = split(input, File.separator);  // Windows "\\", Mac OS "/"
			prefix = String.join(Array.slice(input_split, input_split.length-p, input_split.length), "_");
			if (p != 0) prefix += "_";
			/// Extract file name
			filename = rstrip(file, suffix);  // strip file extension
			/// Combine
			if (remove_num_prefix) sample = prefix + substring(filename, indexOf(filename, "_")+1);  // strip numerical prefix, e.g. 01_t0 -> t0; if not, t0 -> t0, 01_t0 -> 01_t0
			else sample = prefix + filename;
			// Analyze image if this sample has not been completely analyzed (i.e., its output file does not exist)
			if (crop_images) {
				if (!isin(sample + "_crop1_downsize.tif", getFileList(output))) {
					processFile(input, output, file);
				}
			}
			else {
				if (!isin(sample + "_downsize.tif", getFileList(output))) {
					processFile(input, output, file);
				}
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
	// Check if unit is micron
	getPixelSize(unit, pixelWidth, pixelHeight);
	if (unit != "microns") exit("Unit is not microns");
	if (pixelWidth != pixelHeight) exit("pixelWidth != pixelHeight");
	
	// Crop images if needed
	if (crop_images) {
		// Draw crop ROIs (default use the whole image)
		selectWindow(sample); 
		Stack.setChannel(seg_channel_pos);
		run("Select All"); roiManager("add");  // keep whole image as the first ROI
		clear_selection();
		setTool("rectangle");
		waitForUser("Draw one or more rectangular crop regions (include\nsome background) and add to ROI manager [t].\nWhen done, click OK.\nIf not crop, just click OK.");
		clear_selection();
		ncrop = roiManager("count");
		if (ncrop > 1) {  // if crop then more than one ROI, thus delete first ROI (whole image); if not crop then only one ROI, thus use this ROI (whole image)
			roiManager("select", 0); roiManager("delete"); 
			ncrop--; 
			print_and_save_log("Crop " + ncrop + " regions");
		}
		else if (ncrop == 1) print_and_save_log("Not crop");
		clear_selection();
		// Save crop ROIs
		rename_rois();
		roiManager("deselect");
		roiManager("Save", input + File.separator + filename + "_crop_ROIs.zip");
		// Crop and save cropped raw images
		for (c = 1; c <= ncrop; c++) {
			selectWindow(sample);
			roiManager("select", c-1);
			run("Duplicate...", "duplicate");  // crop and duplicate all channels
			saveAs("TIFF", input + File.separator + filename + "_crop" + c + ".tif");
			rename("" + sample + "_crop" + c);
		}
		selectWindow(sample); clear_selection();
		// Clear ROI manager
		roiManager("deselect"); roiManager("delete");
	}
	
	// Process image
	if (!crop_images) ncrop = 1;
	for (c = 1; c <= ncrop; c++) {
		
		// Specify image name
		if (!crop_images) image_name = sample;
		else {image_name = "" + sample + "_crop" + c; print_and_save_log("### Processing crop " + c + " ...");}
		selectWindow(image_name);
		
		// Split channels
		getDimensions(image_width, image_height, nchannel, nslice, nframe);
		if (nchannel == 1) rename("C1-" + image_name);
		else run("Split Channels");
		
		// Cluster segmentation
		// Reason1: use non-cluster region as background region for measuring background intensity for background subtraction during fluorescence quantification
		// Reason2: exclude cluster-center cells which are enriched with overlapping cells, especially for strains with entangled elongated cells
		// Reason3: clear non-cluster region before cell segmentation, because sometimes if there are many debris in non-cluster region, they can significantly affect cellpose cell segmentation
		if (cluster_channel_pos != "") {
			selectWindow("C" + cluster_channel_pos + "-" + image_name);
			if (cluster_channel_pos == seg_channel_pos) run("Duplicate...", " ");
			// Cluster segmentation
			resetMinAndMax(); setOption("ScaleConversions", true); run("8-bit");
			run("Gaussian Blur...", "sigma=0.1 scaled");  //um
			run("Find Edges");
			setAutoThreshold("Li dark");
			setOption("BlackBackground", true);
			run("Convert to Mask");
			run("EDM Binary Operations", "iterations=3 operation=close");  //px
			run("Fill Holes");
			run("Analyze Particles...", "size=" + min_cluster_area + "-Infinity show=Masks in_situ"); // not exclude on edges
			// Remove cluster holes on image edges (artifacts)
			// Algorithm: detect these holes given the fact that they are small and on image edges
			/// Detect edge holes
			cur_image_name = getTitle();
			run("Duplicate...", "title=EdgeHoles");
			run("Invert");
			run("Analyze Particles...", "size=0-" + max_cluster_edgehole_area + " show=Masks in_situ");
			/// Remove non-edge objects
			run("Duplicate...", "title=NonEdge");
			run("Analyze Particles...", "show=Masks exclude in_situ");
			imageCalculator("Subtract", "EdgeHoles", "NonEdge");  // directly modify original image
			/// Fill edge holes
			imageCalculator("Add", cur_image_name, "EdgeHoles");  // directly modify original image
			/// Close temporary images
			selectWindow("EdgeHoles"); close();
			selectWindow("NonEdge"); close();
			// Detect clusters
			run("Analyze Particles...", "size=" + min_cluster_area + "-Infinity show=Masks add in_situ"); // cluster ROIs // not exclude on edges
			clear_selection();
			// Report cluster number
			print_and_save_log("Detect " + roiManager("count") + " clusters");
			// Rename ROIs
			rename_rois();
			// Save cluster ROIs
			roiManager("deselect"); 
			roiManager("save", output + File.separator + image_name + "_cluster_ROIs.zip");
		}
		
		// Downsize and save channel for cell segmentation
		selectWindow("C" + seg_channel_pos + "-" + image_name);
		// Convert to 8-bit
		resetMinAndMax(); setOption("ScaleConversions", true); run("8-bit");
		// Set lookup table to Grays
		run("Grays");
		// Clear non-cluster region
		if (clear_non_cluster) {
			roiManager("Select", Array.getSequence(roiManager("count")));  // cluster ROIs
			roiManager("Combine");
			setBackgroundColor(0, 0, 0);
			run("Clear Outside");
			clear_selection();
		}
		// Save
		saveAs("TIFF", output + File.separator + image_name + "_downsize.tif");  // window name changes to XXX.tif
		
		// Clear ROI manager
		if (roiManager("count") > 0) {roiManager("deselect"); roiManager("delete");}
		
	}

	// Close
	// Close all images
	run("Close All"); 
	// Close ROI manager
	if (!batch_mode_hide) {selectWindow("ROI Manager"); run("Close");}
	// Print complete
	print_and_save_log("Complete processing.");
	print_and_save_log("");	
	
}
