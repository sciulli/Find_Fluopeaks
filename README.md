# Find_Fluopeaks
ImageJ/Fiji tool for automatic detection of fluorescence Ca<sup>2+</sup> transients within x-y-t confocal image stacks. 

##	Description
Automated analysis of intracellular Ca<sup>2+</sup> transients. This ImageJ macro has been developed for the automatic detection of Ca<sup>2+</sup> transients within x-y-t confocal image stacks (design reference formats \*.lif and \*.tif, other formats possibly also compatible) making use of the peak finding plugin tool by Tiago Ferreira (http://fiji.sc/Find_Peaks) included in the BAR collection by the same author (https://imagej.net/BAR). Temporal Ca<sup>2+</sup> transients are evaluated within a user-defined region of interest (using ImageJ selection tools). Fluorescence signals are reported as a ratio (ΔF/F0 %) of the change in fluorescence (ΔF = F-F0) relative to the baseline fluorescence (F0). F0 is automatically detected on the base of the amplitude frequency distribution sampled over a user-selectable number of bins of the fluorescence signal, and hence calculated as the average value (or, optionally, the maximum value) of those signal contributes which are characterized by both high occurence and low amplitude value. The algorithm is designed to detect the Ca<sup>2+</sup> spikes in the time frame between 50 and 250 seconds, which can be modified by the user. Peaks with amplitude lower than standard deviation of (F-F0)/F0\*100 (optionally, 2\*basal signal standard deviation) or value lower than 2\*standard deviation of F0 (condition optionally not considered) are discarded as noise.

##	Requirements
-	ImageJ 1.48d or newer;
-	Find Peaks script by Tiago Ferreira (https://imagej.net/Find_Peaks) which is also included in BAR collection by the same author (https://imagej.net/BAR).

##	Usage
-	Download the macro file;
-	Open the image  stack file with ImageJ/Fiji (tested formats \*.lif and \*.tif, other formats possibly also compatible);
-	Select the region of interest (the algorithm requires a single region to be selected);
-	Open the script and adjust options;
- 	Run the script.

##	Options
- MIN\_ARRAY: analysis starting time point (default value, in sec: 50);
- MAX\_ARRAY: analysis ending time point (default value, in sec: 250);
- N\_BINS: number of histogram bins for baseline signal level evaluation (default value: 20);
- MAX\_DISC_DELTA: minimum numerosity difference needed to discriminate two maxima in the amplitude frequency histogram and therefore contributing to tune the detection accuracy of the algorithm towards the oscillation range of the basal signal (default value: 10);
- HIST\_MAX\_BIN_OFFSET: while considering the amplitude frequency histogram for detecting the fluctuation bin of the basal signal, an offset can be applied in cases where the fluorescence profile presents low-amplitude noisy fluctuations diverse from the actual basal signal (default choice: 0):
	- 0 = consider the basal signal belonging to the lowest bin of the amplitude histogram’s ranked maxima;
	- 1 = consider the basal signal belonging to the first next lowest bin of the amplitude histogram’s ranked maxima;
	- 2 = etc.
- BLINE\_AVG: baseline reference figure (default choice: true):
	- true = average value;
	- false = maximum value.
- MIN\_PEAK\_AMP: peak detection minimum amplitude required (default choice: 1):
	- 1 = F\_norm\_SD; 
	- 2 = bline\_norm\_2SD.
- MIN\_PEAK\_VAL: peak detection minimum amplitude required (default choice: 1): 
	- 1 = bline\_norm\_2SD;
	- 2 = none.

##	Output
Results are shown in ImageJ log and plots windows, including:

###	Basal signal results
-	Basal reference value and number of calculation points;
-	Basal evaluation time interval (bline\_t\_start, bline\_t\_end);
-	Basal oscillation basic statistics: minimum, maximum, mean and standard deviation values (bline\_min, bline\_max, bline\_mean, bline\_SD).

###	Normalized signal results
-	Whole signal basic statistics: minimum, maximum, mean and standard deviation values (F\_norm\_min, F\_norm\_max, F\_norm\_mean, F\_norm_SD);
-	Basal signal basic statistics: minimum, maximum, mean and standard deviation values (bline\_norm\_min, bline\_norm\_max, bline\_norm\_mean, bline\_norm\_SD, bline\_norm\_2SD).

###	Signal peaks results
-	Number of peaks;
-	Frequency;
-	Basic statistics: Mean/Standard Deviation, minimum and maximum values;
-	Coordinates of each peaks in comma-separated form (“t,F\_norm”), easily importable into spreadsheets.

###	Plots
-	Original fluorescence signal (for reference);
-	Fluorescence signal within the selected time interval, detected peaks are marked by red circles.

##	Licensing
This program is free software; you can redistribute and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation (www.gnu.org/licenses/gpl.txt). This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

## Release history note
This (v1.0.1) is the first publicly available release of the script, the very first release dates back to March 2016 when the algorithm started to be employed for analysis of calcium imaging data by Dr. Camilla Fusi, Leibniz Institute for Neurobiology, Magdeburg, Germany.

##	Future developments
Possible further developments include, from usability side, the integration of the algorithm into the ImageJ/Fiji environment in the form of a more conveniently usable plugin tool, provided in particular with UI options panel and more detailed, convenient and versatile results format. On the processing side, the possibility to perform multiple-ROI analysis in a single run may also be added.
