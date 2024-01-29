#@ File (label = "Project directory", style = "directory") proj
#@ String (label = "File suffix", value = ".nd2") suffix
#@ Boolean (label = "Remove numerical prefix of input file names when output (e.g., 01_t0 -> t0)", value = true) remove_num_prefix
#@ String (label = "Channel name(s) (e.g., 'BF eGFP Red')", value = "") channels
#@ String (label = "Measurement(s) (e.g., 'BF:Area,Area raw|eGFP:Mean|Red:Mean,RawIntDen')", value = "") measurements  // see available measurements in ImageJ below
#@ Boolean (label = "Background subtraction for intensity measurements (subtract median background intensity)", value = true) background_subtraction
#@ Boolean (label = "Use batch mode 'hide'", value = true) batch_mode_hide

// Available measurements in ImageJ (see getValue())
// "Area", "Mean", "StdDev", "Mode", "Min", "Max", "X", "Y", "XM", "YM", "Perim.", 
// "BX", "BY", "Width", "Height", "Major", "Minor", "Angle", "Circ.", "Feret", 
// "IntDen", "Median", "Skew", "Kurt", "%Area", "RawIntDen", "Ch", "Slice", "Frame", 
// "FeretX", "FeretY", "FeretAngle", "MinFeret", "AR", "Round", "Solidity", "MinThr", "MaxThr", "Length"

// Script info
pipeline = "cluster_segmentation_measurement";
version = "v2.1";
step = "03";
task = "measurement";
script = String.join(newArray(pipeline, version, step + "_" + task), " ");

// Input/Output
input = proj + File.separator + "00_raw";
input1 = proj + File.separator + "01_auto_segmentation";
input2 = proj + File.separator + "02_manual_correction";
output = proj + File.separator + "03_measurement";
if (!File.exists(input)) exit("Input directory 00_raw does not exist.");
if (!File.exists(input1)) exit("Input directory 01_auto_segmentation does not exist.");
if (!File.exists(input2)) exit("Input directory 02_manual_correction does not exist.");
if (!File.exists(output)) File.makeDirectory(output);

// Print and save log
print_and_save_log("--------------------"); 
print_and_save_log("");
print_and_save_log("Running " + script + " ...");
print_and_save_log("ImageJ version: " + getVersion());
print_and_save_log("OS: " + getInfo("os.name"));
print_and_save_log("Project directory: " + proj);
print_and_save_log("Input directory for raw images: " + input);
print_and_save_log("Input directory for well center ROI (optional): " + input1);
print_and_save_log("Input directory for cluster ROIs: " + input2);
print_and_save_log("Output directory: " + output);
print_and_save_log("File suffix: " + suffix);
print_and_save_log("Remove numerical prefix of input file names when output: " + tf(remove_num_prefix));
print_and_save_log("Channel name(s): " + channels);
print_and_save_log("Measurement(s): " + measurements);
print_and_save_log("Background subtraction for intensity measurements: " + tf(background_subtraction));
print_and_save_log("Use batch mode 'hide': " + tf(batch_mode_hide));
print_and_save_log("Area unit: um^2");
print_and_save_log("");

// Process input parameter values
channels = split(channels, " ");  // e.g. newArray("BF", "eGFP", "Red")
/// Measurements (e.g. "BF:Area,Area raw|eGFP:Mean|Red:Mean,RawIntDen")
measurements_split = split(measurements, "|");
measure_channels = newArray(0);
measure_metrics = newArray(0);
for (i = 0; i < measurements_split.length; i++) {
	x = measurements_split[i];
	x = split(x, ":");
	measure_channels = Array.concat(measure_channels, x[0]);  // e.g. newArray("BF", "eGFP", "Red")
	measure_metrics = Array.concat(measure_metrics, x[1]);  // e.g. newArray("Area,Area raw", "Mean", "Mean,RawIntDen")
}

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
			// Specify input and output file names
			in_ROIs_zip = sample + "_ROIs_auto_manual.zip";
			in_ROIs_eoe_zip = sample + "_ROIs_auto_manual_eoe.zip";
			out_results = sample + "_Results.csv";
			// Analyze image if this sample has not been completely analyzed (i.e., its output file does not exist)
			if (!isin(out_results, getFileList(output))) {
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
	// Check if unit is micron (important for area measurement, etc.)
	getPixelSize(unit, pixelWidth, pixelHeight);
	if (unit != "microns") exit("Unit is not microns");
	if (pixelWidth != pixelHeight) exit("pixelWidth != pixelHeight");
	
	// Load ROIs
	roiManager("Open", input2 + File.separator + in_ROIs_eoe_zip);
	// Print cluster number
	print_and_save_log("Cluster number: " + roiManager("count"));
	
	// Measurement
	// Caveat: manual correction of auto segmentation results may alter some measurements
	run("Clear Results");
	report_background_area = false;  // only need to report background area % once per image
	for (c = 0; c < measure_channels.length; c++) {
		// Set channel
		channel = measure_channels[c];
		if (channels.length > 1) {  // if single channel, then do nothing
			channel_pos = index(channel, channels) + 1;
			Stack.setChannel(channel_pos);
		}
		// Measure
		metrics = measure_metrics[c];
		metrics = split(metrics, ",");
		report_background_intensity = false;  // only need to report background intensity once per channel
		for (m = 0; m < metrics.length; m++) {
			metric = metrics[m];
			colname = "" + channel + "." + metric; 
			// Background subtraction for intensity measurement: measure background area % and median background intensity
			if (background_subtraction && (isin(metric, newArray("Mean", "RawIntDen")))) {
				// Load foreground ROIs
				roiManager("deselect"); roiManager("delete");
				roiManager("Open", input2 + File.separator + in_ROIs_zip);  // object ROIs without excluding on edges
				// Get non-ROI region as background
				roiManager("Select", Array.getSequence(roiManager("count")));
				roiManager("Combine");
				run("Make Inverse");
				// If well center ROI is present, get background within well center
				well_center_present = isin(sample + "_wellcenter.roi", getFileList(input1));
				if (well_center_present) {
					roiManager("Add");  // add "initial" background to ROI manager
					roiManager("Open", input1 + File.separator + sample + "_wellcenter.roi");
					roiManager("select", roiManager("count")-1);
					well_center_area = getValue("Area");
					roiManager("select", newArray(roiManager("count")-2, roiManager("count")-1));
					roiManager("AND");
				}
				// Report background area %
				if (!report_background_area) {
					getDimensions(image_width, image_height, nchannel, nslice, nframe);
					if (well_center_present) background_area_fraction = round(getValue("Area") / well_center_area * 100);
					else background_area_fraction = round(getValue("Area raw") / (image_width*image_height) * 100);
					print_and_save_log("Background area: " + background_area_fraction + "%");
					report_background_area = true;
				}
				// Report median background intensity
				if (!report_background_intensity) {
					median_background_intensity = getValue("Median");
					print_and_save_log("Median background intensity (" + channel + " channel): " + median_background_intensity);
					report_background_intensity = true;
				}
				// Load ROIs with excluding on edges
				roiManager("deselect"); roiManager("delete");
				roiManager("Open", input2 + File.separator + in_ROIs_eoe_zip);
			}
			// Measure
			rn = roiManager("count");
			for (r = 0; r < rn; r++) {
				roiManager("select", r);
				setResult(colname, r, getValue(metric));
				if (background_subtraction && (metric == "Mean")) setResult(colname + "_bgs", r, getValue(metric) - median_background_intensity);
				if (background_subtraction && (metric == "RawIntDen")) setResult(colname + "_bgs", r, getValue(metric) - median_background_intensity * getValue("Area raw"));
			}
			clear_selection();
			// Update results table
			if (!batch_mode_hide) updateResults();
		}
	}
	// Save Results table
	saveAs("Results", output + File.separator + out_results);
	
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
