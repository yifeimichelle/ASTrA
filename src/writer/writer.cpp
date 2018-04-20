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

void write_binned_data(const char *a_filename, int a_numBins, double a_binSize, int a_varDim, const char * const *a_headernames, double **a_vars)
{
  open_file(a_filename);
  // print headers
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

void write_binned_layered_data(const char *a_filename, int a_numBins, double a_binSize, int a_varDim, int a_numLayers, int a_numValues, const char * const *a_headernames, double ****a_vars)
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
      for (unsigned int jPair=0; jPair<a_varDim; jPair++)
  	{
  	  char str[128];
  	  sprintf(str, " %d_0 %d_1", jPair, jPair);
  	  writeString(str);
  	}
      writeString("\n");
      // write data
      for (unsigned int iBin=0; iBin<a_numBins; iBin++)
  	{
  	  char str[128];
  	  sprintf(str, "%f", 1.0*iBin*a_binSize);
  	  writeString(str);
  	  for (unsigned int jPair=0; jPair<a_varDim; jPair++)
  	    {
  	      char str[128];
	      for (int kValue=0; kValue<a_numValues; kValue++)
		{
		  sprintf(str, " %f", a_vars[iLayer][iBin][jPair][kValue]);
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
