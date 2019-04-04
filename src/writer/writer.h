#ifndef _ASTRAWRITER_H_
#define _ASTRAWRITER_H_
#include <iostream>

/// Writes binned data to a file
void write_binned_data(const char *a_filename, int a_numBins, double a_binSize, int a_varDim, const char * const *a_headernames, double **a_vars);
/// Writes binned data with an offset in x-axis to a file
void write_binned_data(const char *a_filename, int a_numBins, double a_binSize, double a_binOffset, int a_varDim, const char * const *a_headernames, double **a_vars);
/// Writes data that is in layers
void write_layered_data(const char *a_filename, int a_numLayers, double *a_layers, int a_varDim, const char * const *a_headernames, double **a_vars);
/// Writes time-series data that is binned and in layers
void write_layered_time_data_by_var(const char *a_filename, int a_numFrames, float a_saveFrameEvery, int a_numLayers, int a_varDim, const char * const *a_headernames, double ***a_vars);
/// Writes time series data
void write_time_data(const char *a_filename, int a_numFrames, float a_saveFrameEvery, int a_varDim, const char * const *a_headernames, double **a_vars);
/// Writes data that is binned and in layers
void write_binned_layered_data(const char *a_filename, int a_numBins, double a_binSize, int a_varDim, int a_numLayers, const char * const *a_headernames, double ***a_vars);
/// Writes data that has bins, layers, and one other dimension (not time)
void write_binned_layered_multival_data(const char *a_filename, int a_numBins, double a_binSize, int a_varDim, int a_numLayers, int a_numValues, const char * const *a_headernames, double ****a_vars);
/// Helper function to write a string
static void writeString(const char *str);
/// Helper function to write a float
static void writeFloat(float val);

#endif
