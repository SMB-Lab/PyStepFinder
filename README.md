# PyStepFinder. 
Author: George Hamilton (ghamil4)
This tool is a Python implementation of a trace data step/change-point finding algorithm. This project was intended as a tool for analysis of fluorescence time-trace data associated with single-molecule FRET experiments and is under active development. The scripts were written in part for participation in the KinSoft challenge (link will be provided upon publication) and the example data included within this repository were generated using the MatLab script provided by the authors of the corresponding publication. Dwell times were calculated from the durations of the steps determined by PyStepFinder, and further analysis was carried out using the dwell times.

At the moment, PyStepFinder is simply a collection of a few functions that allows specifying a test statistic and method of thresholding (difference in RSS, BIC, etc.) and outputs indices of the segment cuts (or a more robust set of outputs if outputstyle='complete').

This early version is provided for use, but potential users should expect drastic changes to the code between versions, including completely restructured code. Additionally, users should be warned that the code is likely inefficient in its calculations, having for the moment been implemented in a naive, brute-force approach.
