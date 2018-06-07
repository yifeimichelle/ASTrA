#include <iostream>
#include <cstring>
#include "writer.h"

static FILE *fp = NULL;

static void open_file(const char *filename)
{
    char full_filename[1024];
    if (strstr(filename, ".out") != NULL)
    {
        strcpy(full_filename, filename);
    }
    else
    {
        sprintf(full_filename, "%s.out", filename);
    }

    fp = fopen(full_filename, "w+");
}

static void close_file(void)
{
  //end_line();
    fclose(fp);
    fp = NULL;
}

// Writes binned data, with multiple variables in a single file
void write_binned_data(const char *a_filename, int a_numBins, double a_binSize, int a_varDim, const char * const *a_headernames, double **a_vars)
{
  open_file(a_filename);
  // print headers
  for (unsigned int iHeader=0; iHeader<a_varDim+1; iHeader++)
    {
      char str[128];
      sprintf(str, " %s", a_headernames[iHeader]);
      writeString(str);
    }
  // print data
  for (unsigned int iBin=0; iBin<a_numBins; iBin++)
    {
      char str[128];
      sprintf(str, "%f", 1.0*iBin*a_binSize);
      writeString(str);
      for (unsigned int jElement=0; jElement<a_varDim; jElement++)
	{
	  char str[128];
	  sprintf(str, " %f", a_vars[iBin][jElement]);
	  writeString(str);
	}
      writeString("\n");
    }
  close_file();
}

// Writes layered data, with multiple variables to a single file
void write_layered_data(const char *a_filename, int a_numLayers, double *a_layers, int a_varDim, const char * const *a_headernames, double **a_vars)
{
  open_file(a_filename);
  // print headers
  for (unsigned int iHeader=0; iHeader<a_varDim+1; iHeader++)
    {
      char str[128];
      sprintf(str, " %s", a_headernames[iHeader]);
      writeString(str);
    }
  writeString("\n");
  // print data
  for (unsigned int iLayer=0; iLayer<a_numLayers; iLayer++)
    {
      char str[128];
      sprintf(str, "%f", a_layers[iLayer]);
      writeString(str);
      for (unsigned int jElement=0; jElement<a_varDim; jElement++)
	{
	  char str[128];
	  sprintf(str, " %f", a_vars[iLayer][jElement]);
	  writeString(str);
	}
      writeString("\n");
    }
  close_file();
}

void write_binned_layered_data(const char *a_filename, int a_numBins, double a_binSize, int a_varDim, int a_numLayers, const char * const *a_headernames, double ***a_vars)
{
  char** filenames = new char * [a_numLayers];
  // compose filenames for each layer
  for (int iLayer=0; iLayer<a_numLayers; iLayer++)
    {
      filenames[iLayer] = new char [1024];
      sprintf(filenames[iLayer], "%s-%d", a_filename, iLayer);
    }

  // Write files for each layer
  for (int iLayer=0; iLayer<a_numLayers; iLayer++)
    {
      open_file(filenames[iLayer]);
      // write headers
      writeString("#r ");
      for (unsigned int jVar=0; jVar<a_varDim; jVar++)
  	{
	  char str[128];
	  sprintf(str, " %d", jVar);
	  writeString(str);
	}
      writeString("\n");
      // write data
      for (unsigned int iBin=0; iBin<a_numBins; iBin++)
  	{
  	  char str[128];
  	  sprintf(str, "%f", 1.0*iBin*a_binSize);
  	  writeString(str);
  	  for (unsigned int jVar=0; jVar<a_varDim; jVar++)
  	    {
	      char str[128];
	      sprintf(str, " %f", a_vars[iLayer][iBin][jVar]);
	      writeString(str);
  	    }
  	  writeString("\n");
  	}
      close_file();
    }
  delete filenames;

}

// Writes time series data, with a different file for each variable
void write_layered_time_data_by_var(const char *a_filename, int a_numFrames, float a_saveFrameEvery, int a_numLayers, int a_varDim, const char * const *a_headernames, double ***a_vars)
{
  char ** filenames = new char * [a_varDim];
  // compose filanemes for each layer
  for (int iVar=0; iVar<a_varDim; iVar++)
    {
      filenames[iVar] = new char [1024];
      sprintf(filenames[iVar], "%s-%d", a_filename, iVar);
    }

  // Write files for each layer
  for (int iVar=0; iVar<a_varDim; iVar++)
    {
      open_file(filenames[iVar]);
      // write headers
      writeString("#t ");
      for (unsigned int jLayer=0; jLayer<a_numLayers; jLayer++)
	{
	  char str[128];
	  sprintf(str, " %d", jLayer);
	  writeString(str);
	}
      writeString("\n");
      // write data
      for (unsigned int iTime=0; iTime<a_numFrames; iTime++)
	{
	  char str[128];
	  sprintf(str, "%f", iTime*a_saveFrameEvery);
	  writeString(str);
	  for (unsigned int jLayer=0; jLayer<a_numLayers; jLayer++)
	    {
	      char str[128];
	      sprintf(str, " %f", a_vars[iTime][jLayer][iVar]);
	      writeString(str);
	    }
	  writeString("\n");
	}
      close_file();
    }
}

// Writes time series data, with multiple variables to a single file
void write_time_data(const char *a_filename, int a_numFrames, float a_saveFrameEvery, int a_varDim, const char * const *a_headernames, double **a_vars)
{
  open_file(a_filename);
  // write headers
  writeString("#t data\n");
  // write data
  for (unsigned int iTime=0; iTime<a_numFrames; iTime++)
    {
      char str[128];
      sprintf(str, "%f", iTime*a_saveFrameEvery);
      writeString(str);
      for (int jVar = 0; jVar<a_varDim; jVar++)
	{
	  sprintf(str, " %f", a_vars[iTime][jVar]);
	  writeString(str);
	}
      writeString("\n");
    }
  close_file();
}



/// Writes binned data with multiple values per bin, to a different file for each layer
void write_binned_layered_multival_data(const char *a_filename, int a_numBins, double a_binSize, int a_varDim, int a_numLayers, int a_numValues, const char * const *a_headernames, double ****a_vars)
{
  char** filenames = new char * [a_numLayers];
  // compose filenames for each layer
  for (int iLayer=0; iLayer<a_numLayers; iLayer++)
    {
      filenames[iLayer] = new char [1024];
      sprintf(filenames[iLayer], "%s-%d", a_filename, iLayer);
    }

  // Write files for each layer
  for (int iLayer=0; iLayer<a_numLayers; iLayer++)
    {
      open_file(filenames[iLayer]);
      // write headers
      writeString("#r ");
      for (unsigned int jVar=0; jVar<a_varDim; jVar++)
  	{
	  for (unsigned int kValue=0; kValue<a_numValues; kValue++)
	    {
	      char str[128];
	      sprintf(str, " %d_%d", jVar, kValue);
	      writeString(str);
	    }
  	}
      writeString("\n");
      // write data
      for (unsigned int iBin=0; iBin<a_numBins; iBin++)
  	{
  	  char str[128];
  	  sprintf(str, "%f", 1.0*iBin*a_binSize);
  	  writeString(str);
  	  for (unsigned int jVar=0; jVar<a_varDim; jVar++)
  	    {
	      for (int kValue=0; kValue<a_numValues; kValue++)
		{
		  char str[128];
		  sprintf(str, " %f", a_vars[iLayer][iBin][jVar][kValue]);
                  writeString(str);
		}
  	    }
  	  writeString("\n");
  	}
      close_file();
    }
  delete filenames;
}

static void writeString(const char *str)
{
    fprintf(fp, "%s",str);
}


static void writeFloat(float val)
{
  char str[128];
  sprintf(str, "%20.12e ", val);
  fprintf(fp, "%s",str);
}
