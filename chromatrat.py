import os;
import re;
import csv;
import string;
import sys;
import math;
import numpy as np;
import matplotlib;
matplotlib.use('Agg');
from pylab import *;

def extractGffFeatures(gffString):
    # GFF v3
    # 1. seqname - Must be a chromosome or scaffold.
    # 2. source - The program that generated this feature.
    # 3. feature - The name of this type of feature. Some examples of standard feature types are "CDS", "start_codon", "stop_codon", and "exon".
    # 4. start - The starting position of the feature in the sequence. The first base is numbered 1.
    # 5. end - The ending position of the feature (inclusive).
    # 6. score - A score between 0 and 1000. If there is no score value, enter ".".
    # 7. strand - Valid entries include '+', '-', or '.' (for don't know/care).
    # 8. frame - If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be '.'.
    # 9. group - All lines with the same group are linked together into a single item.

    res = [];

    if len(gffString) > 1: # check for empty line (i.e. row in GFF file only contains \n)
        try:
            gffSeqName, gffSource, gffFeature, gffStart, gffEnd, gffScr, gffStrand, gffFrame, gffGroup = gffString.split('\t');
        except:
            return []; # exit, since we don't have 9 features

        # handle the gffSeqName
        try:
            chrNo = int(re.findall(r'[0-9]+', gffSeqName)[0]);
        except:
            return [];

        # handle gffStart
        try:
            probeStart = int(gffStart);
        except ValueError:
            return []; # exit, invalid gffStart

        # handle gffEnd
        try:
            probeEnd = int(gffEnd);
        except ValueError:
            return []; # exit, invalid gffEnd

        # handle gffScr
        try:
            probeScr = float(gffScr);
        except ValueError:
            return []; # exit, invalid gffScr

        # if the input string was well-formed, return the extracted features
        res = [chrNo, probeStart, probeEnd, probeScr];
    return res; # return feature set (if any was found in this record)

# Define Constants
tickPlotAdj = 0.5;    # adjust tick marks on x-axis such that they are exactly at the bin boundary (and not in the bin center)

# Initialize Variables
feaNames = [];
feaChrNoList = [];
feaStrandList = [];
feaStartList = [];
feaEndList = [];
feaAttribList = [];
classBoundaryIdxs = [];

gffChrNoList = [];
gffProbeCenterList = [];
gffMatScrList = [];

# get command line arguments
gffFileName = sys.argv[1];
feaFileName = sys.argv[2];
outFileFormat = sys.argv[3];
stepWidth = int(sys.argv[4]);
classBoundaryStr = sys.argv[5];
upStreamStretch = int(sys.argv[6]); # in bp
maxX = int(sys.argv[7]);            # in bp, if -1 maxX <= max(feature length), see below
tickSpacing = int(sys.argv[8]);     # in bp
plotTitle = sys.argv[9];
outFileName = sys.argv[10];

# Extract class boundaries
classBoundaries = map(string.atoi, re.findall(r'\b\d+\b',classBoundaryStr));

# Read transcript information
feaFileReader = csv.reader(open(feaFileName, 'rb'), delimiter = '\t');
for row in feaFileReader:
    feaName, feaChrNo, feaStrand, feaStart, feaEnd, feaAttrib = row;
    feaNames.append(feaName);
    feaChrNoList.append(int(feaChrNo));
    feaStrandList.append(int(feaStrand));
    feaStartList.append(int(feaStart));
    feaEndList.append(int(feaEnd));
    feaAttribList.append(float(feaAttrib));

# Read enrichment scores
gffFileHandle = open(gffFileName, "rU");
for row in gffFileHandle:
    cache = extractGffFeatures(row);
    if len(cache) == 4: # if curr row wasn't well-formed, skip row
        chrNo, probeStart, probeEnd, probeScr = cache;
        gffChrNoList.append(chrNo);
        gffProbeCenterList.append(int(float(probeEnd - probeStart)/2.0 + probeStart)); # calc the center pos of each probe
        gffMatScrList.append(probeScr);

# convert numerical lists into vectors
feaChrNos = np.array(feaChrNoList);
feaStrands = np.array(feaStrandList);
feaStarts = np.array(feaStartList);
feaEnds = np.array(feaEndList);
feaAttribs = np.array(feaAttribList);

gfChrNos = np.array(gffChrNoList);
gfProbeLocs = np.array(gffProbeCenterList);
gfMatScrs = np.array(gffMatScrList);

# preallocate the results array, bins plus 2 extra cols for feature length and transfreq (or other freature characteristic) plus 1 col for sorting
avgFeatureMatScrs = np.empty((len(feaStarts),\
                          (int(math.ceil(float(max(feaEnds - feaStarts))/float(stepWidth)) + upStreamStretch/stepWidth + 3))),\
                         );
avgFeatureMatScrs[:] = np.NAN;

featureCounter = 0;

for chr in range(min(feaChrNos),max(feaChrNos) + 1):
    #get all rows from the feature file that contain information about the i'th chromosome
    currChrIdxs = (feaChrNos == chr);
    currStrands = feaStrands[currChrIdxs];
    currStarts = feaStarts[currChrIdxs];
    currEnds = feaEnds[currChrIdxs];
    currAttribs = feaAttribs[currChrIdxs];

    #get all rows from ChIP-* data sets
    currGfChrIdxs = (gfChrNos == chr);
    currGfProbeLocs = gfProbeLocs[currGfChrIdxs];
    currGfScrs = gfMatScrs[currGfChrIdxs];

    # iterate over all features of the current chromosome
    for j in range(0,len(currStarts)):
        feaStrand = currStrands[j];
        feaStart = currStarts[j];
        feaEnd = currEnds[j];
        feaAttrib = currAttribs[j];

        #calc the corresponding probe positions for the feature start
        if feaStrand == 1:   # +strand
            refFeaSt = 0;
            while feaStart >= currGfProbeLocs[refFeaSt]:
                refFeaSt += 1;
        else:               # -strand
            refFeaSt = 0;
            while feaEnd >= currGfProbeLocs[refFeaSt]:
                refFeaSt += 1;
        #now refFeaSt points to the correct probe according to strand orientation

        #calc the avg probe scores per bin
        if feaStrand == 1:
            feaChunkStart = feaStart;  # bp level pointer
            refChunkStart = refFeaSt;  # probe level pointer
            refChunkEnd = refFeaSt;    # probe level pointer

            # determine values for the upstream bins
            for k in range(1,upStreamStretch/stepWidth + 1):
                feaChunkEnd = feaStart - k * stepWidth;  # bp level pointer

                # find the correct ref for the current feaChunkEnd;  catch potential IndexOutOfBoundEx when upstream region too close to chr boundary
                try:
                    while feaChunkEnd <= currGfProbeLocs[refChunkEnd]:
                        refChunkEnd -= 1;
                except:
                    pass; # don't handle the exception further

                #calc probe mean
                avgFeatureMatScrs[featureCounter, upStreamStretch/stepWidth - k] = np.mean(currGfScrs[refChunkEnd:refChunkStart + 1]);

                #adjust pointers
                refChunkStart = refChunkEnd;

            # determine values for bins in feature body
            for k in range(1, int(math.ceil(float(feaEnd - feaStart)/float(stepWidth)) + 1)):
                feaChunkEnd = feaStart + k*stepWidth;   # bp level pointer

                try:
                    while feaChunkEnd >= currGfProbeLocs[refChunkEnd]:
                        refChunkEnd += 1;
                except:
                    pass; # don't handle the exception further

                # calc probe mean
                avgFeatureMatScrs[featureCounter, k - 1 + upStreamStretch/stepWidth] = np.mean(currGfScrs[refChunkStart:refChunkEnd + 1]);

                # adjust chunk pointers
                refChunkStart = refChunkEnd;

            # write feature length and attrib in last cols of the array
            avgFeatureMatScrs[featureCounter, -1] = feaEnd - feaStart;
            avgFeatureMatScrs[featureCounter, -2] = feaAttrib;
            featureCounter += 1;
        else: # -1 strand
            feaChunkStart = feaEnd;    # bp level pointer
            refChunkStart = refFeaSt;  # probe level pointer
            refChunkEnd = refFeaSt;    # probe level pointer

            # determine values for the upStream region bins
            for k in range(1,upStreamStretch/stepWidth + 1):
                feaChunkEnd = feaEnd + k*stepWidth;  # bp level pointer

                # find the ref for the current feaChunkEnd, catch potential IndexOutOfBoundEx when upstream region too close to chr boundary
                try:
                    while feaChunkEnd >= currGfProbeLocs[refChunkEnd]:
                        refChunkEnd += 1;
                except:
                    pass;

                # calc probe mean
                avgFeatureMatScrs[featureCounter, upStreamStretch/stepWidth - k] = np.mean(currGfScrs[refChunkStart:refChunkEnd+1]);

                # adjust chunk pointers
                refChunkStart = refChunkEnd;

            # determine values for bins in feature body
            for k in range(1, int(math.ceil(float(feaEnd - feaStart)/float(stepWidth)) + 1)):
                feaChunkEnd = feaEnd - k*stepWidth;   # bp level pointer

                try:
                    while feaChunkEnd <= currGfProbeLocs[refChunkEnd]:
                        refChunkEnd -= 1;
                except:
                    pass;

                # calc probe mean
                avgFeatureMatScrs[featureCounter, k-1 + upStreamStretch/stepWidth] = np.mean(currGfScrs[refChunkEnd:refChunkStart+1]);

                # adjust chunk pointers
                refChunkStart = refChunkEnd;

            # write feature length and attrib in last cols of the array
            avgFeatureMatScrs[featureCounter, -1] = feaEnd - feaStart;
            avgFeatureMatScrs[featureCounter, -2] = feaAttrib;
            featureCounter += 1;

# determine idxs of class boundaries and write class value in respective cells
for c in range(len(classBoundaries)-1,-1,-1):   # reverse iterate through the boundaries list
    idxs = np.nonzero(avgFeatureMatScrs[:,-2] <= classBoundaries[c])[0];
    avgFeatureMatScrs[idxs, -3] = c;

# assign names and data format to all cols of the results array and sort by feature attrib
sortByClassCol = 'col' + str(np.size(avgFeatureMatScrs,1)-3);
sortByLengthCol = 'col' + str(np.size(avgFeatureMatScrs,1)-1);

matDtype = {'names':['col%i'%i for i in range(np.size(avgFeatureMatScrs,1))],'formats':(np.size(avgFeatureMatScrs,1))*[np.float]};
avgFeatureMatScrs.dtype = matDtype;
avgFeatureMatScrs.sort(axis = 0, order = [sortByClassCol, sortByLengthCol]);

# generate colormap
BYcmapData = {'red': ((0.0, 0.0, 0.094), (0.5, 0.0, 0.0), (1.0, 0.992, 1.0)), 'green': ((0.0, 0.0,0.658), (0.5, 0.0, 0.0), (1.0, 0.996, 1.0)), 'blue': ((0.0, 1.0, 0.828), (0.5, 0.0, 0.0), (1.0, 0.0, 0.0))};
BYcmap = matplotlib.colors.LinearSegmentedColormap('BYcmap', BYcmapData,9);

lowerHinge = matplotlib.mlab.prctile(gfMatScrs, 25.0);
upperHinge = matplotlib.mlab.prctile(gfMatScrs, 75.0);
cmLim = (upperHinge - lowerHinge) * 1.5 + upperHinge;   # upper inner fence

# remove custom array dtype used for sorting and plot avg scrs (without the last 3 cols)
imgplot = imshow(avgFeatureMatScrs.view(np.float)[:,0:-3], aspect='auto', interpolation='nearest', origin='lower', vmin = -cmLim, vmax = cmLim);
axvline(x = int(float(upStreamStretch)/float(stepWidth)) - tickPlotAdj, color = 'w', linewidth = 1);

# if maxX wasn't set by hand, set it to max feature length and set x-coord accordingly
if maxX == -1:
    maxX = max(feaEnds - feaStarts);
    imgplot.axes.set_xbound(-tickPlotAdj, int(float(maxX)/float(stepWidth)) + int(float(upStreamStretch)/float(stepWidth)) - tickPlotAdj);
else:
    imgplot.axes.set_xbound(-tickPlotAdj, int(float(maxX)/float(stepWidth)) - tickPlotAdj);
    
xTicks = [k * int(float(tickSpacing)/float(stepWidth)) + int(float(upStreamStretch)/float(stepWidth)) - tickPlotAdj for k in range(0, int(math.ceil(float(maxX)/float(tickSpacing)))+1)];
xTicks.insert(-1,-tickPlotAdj);
xTicksLabel = [k * tickSpacing for k in range(0, int(math.ceil(float(maxX)/float(tickSpacing)))+1)];
xTicksLabel.insert(-1, str(-upStreamStretch));
imgplot.axes.set_xticks(xTicks);
imgplot.axes.set_xticklabels(xTicksLabel);

imgplot.axes.set_title(plotTitle);
imgplot.axes.set_xlabel('Distance from feature start (bp)');
imgplot.axes.set_ylabel('Feature class (low to high)');
imgplot.axes.set_yticks([]);

imgplot.set_cmap(BYcmap);
colorbar();
show();

savefig('chromatra_t_tmp_plot', format=outFileFormat);
data = file('chromatra_t_tmp_plot', 'rb').read();
fp = open(outFileName, 'wb');
fp.write(data);
fp.close();
os.remove('chromatra_t_tmp_plot');
