
/* Find_Fluopeaks.ijm
 *
 * Find Fluopeaks, ImageJ/Fiji tool for automatic detection of fluorescence Ca2+ transients within x-y-t confocal image stacks.
 * Requires ImageJ 1.48d or newer.
 * Requires Find Peaks script by Tiago Ferreira (https://imagej.net/Find_Peaks) which is also included in BAR collection by the same author (https://imagej.net/BAR).
 * See https://github.com/sciulli/Find_Fluopeaks for more details.
 * Stefano Ciulli, v1.0.1 2019.11.25
 * <sciulli.dev (at) gmail.com.nospam>
 * 
 * LICENSING:
 * This program is free software; 
 * you can redistribute and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation (www.gnu.org/licenses/gpl.txt). 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
 * See the GNU General Public License for more details. 
 */

// Set debugging flag
DEBUG=false;
CLOSE_ALL_PLOTS=true;
CLOSE_PEAK_PLOT=false;

// Processing variables
MIN_ARRAY = 50; // analysis starting time point
MAX_ARRAY = 250; // analysis ending time point
N_BINS = 20; // number of histogram bins for baseline signal level evaluation
MAX_DISC_DELTA = 10; // minimum numerosity difference needed to discriminate two maxima in the amplitude frequency histogram and therefore contributing to tune the detection accuracy of the algorithm towards the oscillation range of the basal signal
HIST_MAX_BIN_OFFSET = 0; //offset possibly applied in cases where the fluorescence profile presents low-amplitude noisy fluctuations diverse from the actual basal signal
// 0 = consider the basal signal belonging to the lowest bin of the amplitude histogram’s ranked maxima;
// 1 = consider the basal signal belonging to the first next lowest bin of the amplitude histogram’s ranked maxima;
// 2 = etc.;
BLINE_AVG = true; // baseline reference figure, true -> average value, false -> maximum value;
MIN_PEAK_AMP = 1; // peak detection minimum amplitude required, 1 -> F_norm_SD; 2 -> bline_norm_2SD
MIN_PEAK_VAL = 1; // peak detection minimum value required, 1 -> bline_norm_2SD; 2 -> none

imgID = getImageID(); //keep the ID of the working image
getSelectionBounds(roi_x, roi_y, roi_width, roi_height); //get user defined roi position and size

// Show original fluorescence profile
run("Plot Z-axis Profile");
Plot.getValues(X_, Y_); // values are stored into "Array"
X = Array.slice(X_,MIN_ARRAY,MAX_ARRAY);
Y = Array.slice(Y_,MIN_ARRAY,MAX_ARRAY);
//if (CLOSE_ALL_PLOTS) {
//	close(); // close also original signal plot window (image type window)
//}	
closeWin(); // close Results window

// Show fluorescence profile within the selected time period (MIN_ARRAY, MAX_ARRAY)
plotTitle = "Original signal truncated between " + MIN_ARRAY + " and " + MAX_ARRAY + " sec";
Plot.create(plotTitle, "X [sec]", "F" , X, Y);
Plot.show();
if (DEBUG) {
	for (i=0; i<X.length; i++) {
		print("X[" + i + "]: " + X[i], "Y[" + i + "]: " + Y[i]);
	}
}
if (CLOSE_ALL_PLOTS){
	close(); // close plot window (image type window)
}
// Get statistics from data array
Array.getStatistics(X, X_min, X_max, X_mean, X_stdDev);
Array.getStatistics(Y, Y_min, Y_max, Y_mean, Y_stdDev);
// Arrange bins for building profile histogram 
Nbins = N_BINS;
bin_size = floor(X_max/Nbins) + 1;
bin_range = floor(Y_max/Nbins) + 1;
bin_min = newArray(Nbins);
bin_max = newArray(Nbins);
totcount = 0;
for (b=0; b<Nbins; b++) {
	bin_min[b] = b*bin_range;
	bin_max[b] = (b+1)*bin_range;
}
// Build profile amplitude's frequency histogram 
counts = newArray(Nbins);
for (i=0; i<Y.length; i++) {
	for (b=0; b<Nbins; b++) {
  		if (bin_min[b] <= Y[i] && Y[i] < bin_max[b]) {  
  			counts[b]++;
  			totcount++;
  		}
	}
}
// Show frequency histogram values for debugging purposes
if (DEBUG) {
	for (b=0; b<Nbins; b++) {
		print(b + " " + bin_min[b] + " " + bin_max[b] + " " + counts[b]);
	}
	print(Nbins + " " + Y_max + " " + bin_range + " " + totcount);
}

// Build arrays holding maxima values and frequences of the profile amplitude's frequency histogram 
maxLocs = Array.findMaxima(counts, MAX_DISC_DELTA);
hist_max_x = newArray(maxLocs.length);
hist_max_y = newArray(maxLocs.length);
for (jj = 0; jj < maxLocs.length; jj++){
	hist_max_x[jj] = maxLocs[jj];
	hist_max_y[jj] = counts[hist_max_x[jj]];
	if (DEBUG) {
		print("hist_max_x = ", hist_max_x[jj], " hist_max_y = ", hist_max_y[jj]);
	}
}

// Build an array holding the rank position indexes of maxima values of the profile amplitude's frequency histogram, starting with the index of the smallest VALUE
hist_max_x_rank = Array.rankPositions(hist_max_x); 
// Get the index of the lowest(*) bin of amplitude frequency histogram maxima, considered to be the bin which the basal signal belongs to;
// (*) an offset, represented by HIST_MAX_BIN_OFFSET option value, can be applied in cases where the fluorescence profile presents low-amplitude noisy fluctuations diverse from the actual basal signal
hist_max_x_minrank = hist_max_x_rank[HIST_MAX_BIN_OFFSET]; 

// Basal signal bin minimum and maximum values
bline_min = bin_min[hist_max_x[hist_max_x_minrank]];
bline_max = bin_max[hist_max_x[hist_max_x_minrank]];
if (DEBUG) {
	print("bline_min = ", bline_min, " bline_max = ", bline_max);
}

// Scanning all profile points within the baseline bin: [bline_min, bline_max], corresponding coordinate arrays are built
baseline = 0;
NblinePoints = -1; //(incremental) number of baseline calculation points
NavgPoints = hist_max_y[hist_max_x_minrank]; //number of baseline calculation points
X_bline = newArray(NavgPoints);
Y_bline = newArray(NavgPoints);
indexes_bline = newArray(NavgPoints);
for (i=0; i<Y.length; i++) {
	if (Y[i] >= bline_min && Y[i] < bline_max) {
		NblinePoints++;
		X_bline[NblinePoints] = X[i];
		Y_bline[NblinePoints] = Y[i];
		indexes_bline[NblinePoints] = i;
	}
}

// Basal signal statistics
Array.getStatistics(X_bline, X_bline_min, X_bline_max, X_bline_mean, X_bline_stdDev);
Array.getStatistics(Y_bline, Y_bline_min, Y_bline_max, Y_bline_mean, Y_bline_stdDev);

// Baseline value assignment in accordance to user defined options
if (BLINE_AVG){
	baseline = Y_bline_mean; // baseline average value
}
else{
	baseline = Y_bline_max; // baseline maximum value
}

// Baseline detection results output
print("**** BASAL SIGNAL RESULTS ****");
print("baseline:", baseline, "calculated on", NavgPoints, "points.");
print("bline_t_start:", X_bline_min, "\nbline_t_end:", X_bline_max);
print("bline_min:", Y_bline_min, "\nbline_max:", Y_bline_max, "\nbline_mean:", Y_bline_mean, "\nbline_SD:", Y_bline_stdDev);
print(""); // emptyline

// Calculate normalized signal (F-F0)/F0 * 100
Y_norm = newArray(Y.length);
for (i=0; i<Y.length; i++) {
	Y_norm[i] = ((Y[i] - baseline)/baseline)*100;
	if (DEBUG) {
		print("X = ", X[i], " Y_norm = ", Y_norm[i]);	
	}
}
// Nomalized signal statistics
Array.getStatistics(Y_norm, Y_norm_min, Y_norm_max, Y_norm_mean, Y_norm_stdDev);

// Plot normalized sgnal 
Plot.create("F_norm [(F-F0)/F0 %]", "t [sec]", "F_norm [(F-F0)/F0 %]", X, Y_norm);
Plot.show();

// Isolate baseline normalized data points 
Y_bline_norm = newArray(NavgPoints);
for (b=0; b<X_bline.length; b++){
	Y_bline_norm[b] = Y_norm[indexes_bline[b]];
}
// basal signal normalized data points statistics
Array.getStatistics(Y_bline_norm, Y_bline_norm_min, Y_bline_norm_max, Y_bline_norm_mean, Y_bline_norm_stdDev);
Y_bline_norm_2SD = Y_bline_norm_stdDev*2;

// Normalized signal results output
print("**** NORMALIZED SIGNAL RESULTS ****");
print(">> normalized fluorescence whole signal basic statistics");
print("F_norm_min:", Y_norm_min, "\nF_norm_max:", Y_norm_max, "\nF_norm_mean:", Y_norm_mean, "\nF_norm_SD:", Y_norm_stdDev);
print(">> normalized basal signal basic statistics");
print("bline_norm_min:", Y_bline_norm_min, "\nbline_norm_max:", Y_bline_norm_max, "\nbline_norm_mean:", Y_bline_norm_mean, "\nbline_norm_SD:", Y_bline_norm_stdDev, "\nbline_norm_2SD:", Y_bline_norm_2SD);
print(""); // emptyline

// Set Find Peaks command line options in accordance to the user defined options  
if (MIN_PEAK_AMP == 1){
	MIN_PEAK_AMPLI = "Y_norm_stdDev"; // minimum peak amplitude = normalized fluorescence standard deviation
}
else if (MIN_PEAK_AMP == 1){
	MIN_PEAK_AMPLI = "Y_bline_norm_2SD"; // minimum peak amplitude = 2 * normalized baseline standard deviation
}
if (MIN_PEAK_VAL == 1) {
	MIN_PEAK_VALUE = "Y_bline_norm_2SD"; // minimum peak value = 2 * normalized baseline standard deviation
}
else if (MIN_PEAK_VAL == 2) { // minimum peak value = none
	MIN_PEAK_VALUE = "[]";
}

// Run Find Peaks plugin by Tiago Ferreira (https://imagej.net/Find_Peaks)
run("Find Peaks", "min._peak_amplitude="+MIN_PEAK_AMPLI+" min._peak_distance=0 min._value="+MIN_PEAK_VALUE+" max._value=[] list");
selectWindow("Plot Values");
NresRows = getValue("results.count");
peaks_X = newArray(NresRows);
peaks_Y = newArray(NresRows);
peaksmin_X = newArray(NresRows);
peaksmin_Y = newArray(NresRows);
// peaks_HM = newArray(NresRows); //half-maximum, not used
Npeaks = 0;
Npeaksmin = 0;
// Get detected peaks coordinates
for (row=0; row < NresRows; row++) {
	if (!isNaN(getResult("Y1",row))){ // peaks top
		peaks_X[row] = getResult("X1",row);
		peaks_Y[row] = getResult("Y1",row);
		// peaks_HM[row] = peaks_Y[row]/2; //half-maximum, not used
		Npeaks++;
	}
	if (!isNaN(getResult("Y2",row))){ // peaks bottom
		peaksmin_X[row] = getResult("X2",row);
		peaksmin_Y[row] = getResult("Y2",row);
		Npeaksmin++;
	}
}
// Turn row arrays into column arrays
peaks_X = Array.slice(peaks_X,0,Npeaks);
peaks_Y = Array.slice(peaks_Y,0,Npeaks);
//peaks_HM = Array.slice(peaks_HM,0,Npeaks); // peaks half-maximum, not used

// Close some unused windows
selectWindow("Plot Values");
run("Close");
if (CLOSE_ALL_PLOTS){
	selectWindow("F_norm [(F-F0)/F0 %]");
	close();
	if (CLOSE_PEAK_PLOT){
		selectWindow("Peaks in F_norm [(F-F0)/F0 %]");
		close();
	}
}

// Peaks data points statistics
Array.getStatistics(peaks_X, peaks_X_min, peaks_X_max, peaks_X_mean, peaks_X_stdDev);
Array.getStatistics(peaks_Y, peaks_Y_min, peaks_Y_max, peaks_Y_mean, peaks_Y_stdDev);
frequency = peaks_X.length / (peaks_X_max-peaks_X_min);

// Peaks details output
print("**** PEAKS RESULTS ****");
print("N. peaks: " + Npeaks);
print("Frequency: " + frequency + " [pks/sec] over " + (peaks_X_max-peaks_X_min) + " sec.");
print("Mean/SD amplitude:",peaks_Y_mean,"/",peaks_Y_stdDev,"\nMin amplitude:",peaks_Y_min,"\nMax amplitude:",peaks_Y_max);
print("Peaks(t,F_norm):");
print("t,F_norm");
for (p=0; p < Npeaks; p++) {
//print(peaks_X[p]+","+peaks_Y[p]+","+peaks_HM[p]); // include half-maximum, not used
print(peaks_X[p]+","+peaks_Y[p]);
}

exit()

function closeWin() {
	winList = getList("window.titles");
	for (w=0; w<winList.length; w++){
		//print(winList[w]); // test purpose
		ROImanagerTitle = "ROI Manager";
		LogTitle = "Log";
		if ((winList[w] != ROImanagerTitle)&&(winList[w] != LogTitle)){
			selectWindow(winList[w]);
			run("Close");
		}
	}
}

/*
 * LICENSING:
 * This program is free software; 
 * you can redistribute and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation (www.gnu.org/licenses/gpl.txt). 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
 * See the GNU General Public License for more details. 
 */

