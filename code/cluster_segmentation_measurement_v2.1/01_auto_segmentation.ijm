#@ File (label = "Project directory", style = "directory") proj
#@ String (label = "File suffix", value = ".nd2") suffix
#@ Boolean (label = "Remove numerical prefix of input file names when output (e.g., 01_t0 -> t0)", value = true) remove_num_prefix
#@ Integer (label = "Position of BF channel for cluster segmentation", value = 1) BF_pos
#@ Float (label = "Min cluster area (um^2)", value = 90) min_cluster_area
#@ Boolean (label = "Whole-well image", value = false) whole_well_image
#@ Boolean (label = "If auto well center (w/wo touching clusters) detection fails: skip file (uncheck), manually correct (check)", value = false) manual_correct_auto_well_center
#@ Float (label = "Shrink well center radius (um) after auto well center detection (e.g., 400, 0)", value = 400) shrink_well_center_radius
#@ Boolean (label = "Include clusters touching left (but not right) border of well center, otherwise exclude both", value = true) include_wellcenter_left
#@ Boolean (label = "Use batch mode 'hide'", value = true) batch_mode_hide

// Script info
pipeline = "cluster_segmentation_measurement";
version = "v2.1";
step = "01";
task = "auto_segmentation";
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
print_and_save_log("Position of BF channel for cluster segmentation: " + BF_pos);
print_and_save_log("Min cluster area (um^2): " + min_cluster_area);
print_and_save_log("Whole-well image: " + tf(whole_well_image));
print_and_save_log("If auto well center (w/wo touching clusters) detection fails: skip file (uncheck), manually correct (check): " + tf(manual_correct_auto_well_center));
print_and_save_log("Shrink well center radius (um) after auto well center detection: " + shrink_well_center_radius);
print_and_save_log("Include clusters touching left (but not right) border of well center, otherwise exclude both: " + tf(include_wellcenter_left));
print_and_save_log("Use batch mode 'hide': " + tf(batch_mode_hide));
print_and_save_log("");

// Parameters
alt_radius = 50;  //px  // auto local threshold (alt)
alt_par1 = 0.15;  // default 0.25 in ImageJ's auto local threshold command with Phansalkar method 
alt_par2 = 0.5;  // default 0.5 in ImageJ's auto local threshold command with Phansalkar method 
edm_close_cycles = 4;
gaussian_filter_radius = 7;  //um
find_maxima_prominence = 5;

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
			if (remove_num_prefix) filename = substring(filename, indexOf(filename, "_")+1);  // strip numerical prefix, e.g. 01_t0 -> t0; if not, t0 -> t0, 01_t0 -> 01_t0
			/// Combine
			sample = prefix + filename;
			// Analyze image if this sample has not been completely analyzed (i.e., its final output file does not exist)
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
	// Check if unit is micron (important for area thresholds, etc.)
	getPixelSize(unit, pixelWidth, pixelHeight);
	if (unit != "microns") exit("Unit is not microns");
	if (pixelWidth != pixelHeight) exit("pixelWidth != pixelHeight");
	
	// Split channels
	getDimensions(image_width, image_height, nchannel, nslice, nframe);
	if (nchannel == 1) rename("C1-" + file);
	else run("Split Channels");
	// Close channels other than BF
	images = getList("image.titles");
	for (i = 0; i < images.length; i++) {
		selectWindow(images[i]);
		if (images[i] != "C" + BF_pos + "-" + file) close();
	}
	
	// Cluster segmentation
	// Convert to 8-bit
	resetMinAndMax(); setOption("ScaleConversions", true); run("8-bit");
	// Save downsized image for manual correction
	saveAs("TIFF", output + File.separator + sample + "_downsize.tif");  // window name changes to XXX.tif
	rename("Raw");
	// Segmentation
	run("Duplicate...", "title=Mask");
	run("Auto Local Threshold", "method=Phansalkar radius=" + alt_radius + " parameter_1=" + alt_par1 + " parameter_2=" + alt_par2);
	if (whole_well_image) {
		// Detect bright well center (excluding dark well edge region if present)
		run("Invert");
		run("Analyze Particles...", "size=0-Infinity show=Masks exclude in_situ");  // exclude on edges (necessary for some multi-well plates)
		run("Analyze Particles...", "size=" + (image_width * image_height * 0.3) + "-Infinity pixel exclude include add");
		if (!detect_well_center("")) return;  // stop processing this file and move to the next file
		if (shrink_well_center_radius != 0) {  // if == 0 then not shrink well center radius (optional if no dark well edge), essentially detect well border
			roiManager("Select", 0);
			run("Enlarge...", "enlarge=-" + shrink_well_center_radius);  // shrink well center by a little to exclude regions that are potentially dark or may affect cluster segmentation or fluorescence quantification
			roiManager("Update");
		}
		roiManager("Select", 0);
		roiManager("rename", "wellcenter");
		roiManager("Save", output + File.separator + sample + "_wellcenter.roi");
		clear_selection();
		run("Invert");
		// Clear outside of well center but retain clusters touching border of well center
		if (shrink_well_center_radius != 0) {  // if == 0 then not shrink well center radius (optional if no dark well edge), essentially clear outside of well border
			run("Duplicate...", "title=ClearOutside");
			roiManager("Select", 0);
			setForegroundColor(255, 255, 255);
			run("Fill", "slice");
			clear_selection();
			roiManager("deselect"); roiManager("delete");
			run("Analyze Particles...", "size=" + (image_width * image_height * 0.3) + "-Infinity pixel exclude add");
			if (!detect_well_center("(with touching clusters) ")) return;  // stop processing this file and move to the next file
			clear_selection();
			close();
		}
		roiManager("Select", 0);
		setBackgroundColor(0, 0, 0);
		run("Clear Outside");
		clear_selection();
		roiManager("deselect"); roiManager("delete");
	}
	run("EDM Binary Operations", "iterations=" + edm_close_cycles + " operation=close");
	run("Analyze Particles...", "size=" + min_cluster_area + "-Infinity show=Masks add in_situ");  // add to ROI manager for "Add back small objects whose maxima were not detected" below // filtering based on min cluster area is fine even before Fill Holes because those tiny objects do not need Fill Holes
	clear_selection();
	// Split touching clusters
	/// Generate seed image
	selectWindow("Raw");
	run("Duplicate...", "title=Seed");
	run("Gaussian Blur...", "sigma=" + gaussian_filter_radius + " scaled");
	run("Find Maxima...", "prominence=" + find_maxima_prominence + " strict light output=[Single Points]");  // generate "* Maxima" image  // maxima of some small objects (tiny propagules and debris) may not be detected, add back later
	selectWindow("Seed"); close();
	/// Watershed of masks based on seeds and raw image
	/// Note: if no seed but mask then mask is not shown, if seed but no mask then mask is also not shown
	run("Marker-controlled Watershed", "input=Raw marker=[Seed Maxima] mask=[Mask] compactness=0 binary calculate use");  // generate 32-bit image (labeled masks), but not all values between 1-max have corresponding pixel values or masks
	rename("Clusters");
	setOption("ScaleConversions", false); run("8-bit"); setOption("ScaleConversions", true);  // all pixel values above 255 will be 255, otherwise keep original value
	setThreshold(1, 255); setOption("BlackBackground", true); run("Convert to Mask");
	/// Add back small objects whose maxima were not detected
	for (r = 0; r < roiManager("count"); r++) {
	    roiManager("select", r);
	    if (getValue("RawIntDen") == 0) {setForegroundColor(255, 255, 255); run("Fill", "slice");}
	}
	clear_selection();
	roiManager("deselect"); roiManager("delete"); 
	/// Close temporary images
	selectWindow("Raw"); close();
	selectWindow("Seed Maxima"); close();
	selectWindow("Mask"); close();
	// Fill holes
	selectWindow("Clusters");
	run("Fill Holes");  // fill holes after splitting touching clusters, to handle situations where multiple clusters surround a background area
	if (whole_well_image) {
		// For clusters touching border of well center, either exclude all of them or only keep those touching left side of border
		// If not shrink well center radius before (optional if no dark well edge), then exclude all of them
		/// Load well center ROI
		roiManager("Open", output + File.separator + sample + "_wellcenter.roi");
		wc = roiManager("count") - 1; // ROI index of well center ROI
		roiManager("select", wc);
		roiManager("rename", "wellcenter");
		/// If include clusters touching the left border of well center
		if (include_wellcenter_left) {
			// Split bright well center ROI into left and right halves
			/// Get x coordinate of the center of bright well center 
			roiManager("select", wc);
			getPixelSize(unit, pixelWidth, pixelHeight);
			wc_center_x = round(getValue("X") / pixelWidth);
			/// Save middle line ROI that splits well center into left/right halves
			getDimensions(image_width, image_height, nchannel, nslice, nframe);  //px
			makeLine(wc_center_x, 0, wc_center_x, image_height);
			roiManager("Add");
			midline = roiManager("count") - 1;  // ROI index
			roiManager("Select", midline);
			roiManager("rename", "wellcenter_midline");
			roiManager("Save", output + File.separator + sample + "_wellcenter_midline.roi");
			/// Create left rectangle ROI
			makeRectangle(0, 0, wc_center_x, image_height);
			roiManager("Add");
			rect = roiManager("count") - 1;  // ROI index
			/// Get left half of well center
			roiManager("Select", newArray(wc, rect));
			roiManager("AND");
			roiManager("Add");
			wc_left = roiManager("count") - 1;  // ROI index
			/// Get right half of well center
			roiManager("Select", newArray(wc, wc_left)); 
			roiManager("XOR");
			roiManager("Add");
			wc_right = roiManager("count") - 1;  // ROI index
			clear_selection();
			// Expand left well center ROI to include cluster ROIs that touch the border of left well center ROI
			run("Duplicate...", "title=IncludeLeft");
			roiManager("select", wc_left);
			wc_left_area = getValue("Area");
			setForegroundColor(255, 255, 255);
			run("Fill", "slice");
			clear_selection();
			run("Analyze Particles...", "size=" + (wc_left_area-1) + "-Infinity show=Nothing add in_situ"); 
			clear_selection();
			wc_left_expand = roiManager("count") - 1;  // ROI index
			if (wc_left_expand - wc_left != 2) exit("ROI index: wc_left_expand - wc_left != 2");
			selectWindow("IncludeLeft"); close();
			// Combine expanded left well center and right well center
			roiManager("select", wc_left_expand);
			run("Enlarge...", "enlarge=1 pixel");  // all non-touching clusters should be away by at least 1px
			roiManager("Add");
			wc_left_expand_plus1px = roiManager("count") - 1;  // ROI index
			roiManager("Select", newArray(wc_left_expand_plus1px, rect));
			roiManager("AND");
			roiManager("Add");
			wc_left_expand_plus1px_trim = roiManager("count") - 1;  // ROI index
			roiManager("Select", newArray(wc_left_expand_plus1px_trim, wc_right));
			roiManager("Combine");
			roiManager("Add");
			wc_expand = roiManager("count") - 1;  // ROI index
			roiManager("Select", wc_expand);
		}
		/// Exclude outside clusters
		run("Make Inverse");
		setForegroundColor(255, 255, 255);
		run("Fill", "slice");
		clear_selection();
		run("Analyze Particles...", "size=0-Infinity show=Masks exclude in_situ");  // exclude on edges
		/// Delete all intermediate ROIs
		roiManager("deselect"); roiManager("delete"); 
	}
	run("Analyze Particles...", "size=" + min_cluster_area + "-Infinity show=Masks in_situ add");  // not exclude on edges here but do this later after manual correction
	// Report cluster number
	print_and_save_log("Cluster number after auto segmentation: " + roiManager("count"));
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

// Check if single well center is detected
// Detected well center (ROI number in ROI manager) should be 1 but may be 0 or more than 1, in the latter cases, perform manual correction
// well_center_type = "", "(with touching clusters) "
function detect_well_center(well_center_type) { 
	if (!manual_correct_auto_well_center) {
		if (roiManager("count") == 1) is_well_center_detected = true;
		else is_well_center_detected = false;
	}
	else {
		while (roiManager("count") != 1) {  
			clear_selection();
			setBatchMode("show");
			waitForUser("" + roiManager("count") + " well centers " + well_center_type + "are detected.\n" +
						"Manually correct (e.g., select none, invert, manually separate well center and others using pencil tool, select none, invert).");
			setBatchMode("hide");
			clear_selection();
			print_and_save_log("Manually correct to detect single well center " + well_center_type);
			if (roiManager("count") >= 1) {roiManager("deselect"); roiManager("delete");}
			run("Analyze Particles...", "size=" + (image_width * image_height * 0.3) + "-Infinity pixel exclude include add");
		}
		if (roiManager("count") == 1) is_well_center_detected = true;
		else is_well_center_detected = false;
	}
	if (!is_well_center_detected) {
		print_and_save_log("Single well center " + well_center_type + "is not detected. Skip this file.");
		print_and_save_log("");
	}
	return is_well_center_detected;
}
