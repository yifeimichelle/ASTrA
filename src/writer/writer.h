#ifndef _ASTRAWRITER_H_
#define _ASTRAWRITER_H_
#include <iostream>

/// Writes binned data to a file
void write_binned_data(const char *a_filename, int a_numBins, double a_binSize, int a_varDim, const char * const *a_headernames, double **a_vars);
void write_binned_layered_data(const char *a_filename, int a_numBins, double a_binSize, int a_varDim, int a_numLayers, int a_numValues, const char * const *a_headernames, double ****a_vars);
static void writeString(const char *str);
static void writeFloat(float val);

#endif
