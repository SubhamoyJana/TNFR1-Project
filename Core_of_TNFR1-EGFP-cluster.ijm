//This code segments TNFR1-EGFP clusters into core and rim regions based on intensity thresholds. 
//It categorizes pixels with intensities greater than 70% of the maximum intensity as core and 
//those between 10% and 70% as rim. You should open Z-stack images of individual TNFR1-EGFP clusters
//with total intensity and anisotropy, duplicate the total intensity image, correct the background, 
//and then run the code.
run("Clear Results");
activesite = getTitle();
selectWindow(activesite);
run("32-bit");
w = getWidth(); h = getHeight();
n = nSlices;
run("Set Measurements...", "area mean min redirect=None decimal=3");
for (i = 1; i < n+1; i++) {
	selectWindow(activesite);
	setSlice(i);
	//roiManager("select",i-1);
	run("Measure");
	MaxInt = getResult("Max", i-1);
	MinInt = getResult("Min", i-1);
	run("Divide...", "value="+MaxInt+" slice");
	}
selectWindow(activesite);
//setThreshold(0.1,0.7);  // for rim
setThreshold(0.7, 1); // for core
run("NaN Background", "stack");
// after getting the threshold image take the ROI by doWand tool. That will be the core region
	
