#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix", value = "_Results.csv") suffix

processFolder(input);

function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i]);
	}
}

function processFile(input, output, file) {
	open(input + File.separator + file);
	Table.rename(file, "Results");
	Plot.create("Plot of Results", "Area", "cell_top5_mean_area");
	Plot.add("Circle", Table.getColumn("Area", "Results"), Table.getColumn("cell_top5_mean_area", "Results"));
	Plot.setStyle(0, "blue,#a0a0ff,1.0,Circle");
	Plot.setLimits(0, 21000, 0, 26);
	Plot.show();
	saveAs("PNG", output + File.separator + rstrip(file, suffix) + "_PlotLimits.png");
	run("Clear Results");
	run("Close All");
}

function rstrip(string, suffix) {
	return substring(string, 0, lengthOf(string) - lengthOf(suffix));	
}
