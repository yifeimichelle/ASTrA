#include "frame.h"
#include "system.h"
#include "atom.h"
#include "radialDistribution.h"
#include <vector>
#include <iostream>
#include <assert.h>

using namespace std;

static FILE *fp = NULL;

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

RDF::RDF()
{
};

RDF::RDF(System& a_system)
{
  m_system = a_system;
  m_numLayers = m_system.getNumLayers();
  m_rdf.resize(m_system.getNumPairs());
  m_rdfLayer.resize(m_numLayers);
  m_rdfLayerClosest.resize(m_numLayers);

  m_pairCounter.resize(m_system.getNumPairs());
  m_maxDist = 14.0;
  m_numBins = 500;
  for (int i=0; i<m_system.getNumPairs(); i++)
    {
      m_rdf[i].resize(m_numBins);
    }
  for (int i=0; i<m_numLayers; i++)
    {
      m_rdfLayer[i].resize(m_system.getNumPairs());
      m_rdfLayerClosest[i].resize(m_system.getNumPairs());
      for (int j=0; j<m_system.getNumPairs(); j++)
	{
	  m_rdfLayer[i][j].resize(m_numBins);
	  m_rdfLayerClosest[i][j].resize(2);
	  for (int k=0; k<2; k++)
	    {
	      m_rdfLayerClosest[i][j][k].resize(m_numBins);
	    }
	}
    }
  m_binSize = m_maxDist / m_numBins;
};

void RDF::sample(Frame& a_frame)
{
  double pairDistance;
  double minDistance;
  for (int pairIdx=0; pairIdx<m_system.getNumPairs(); pairIdx++) {
    pair<unsigned int, unsigned int > tpair = m_system.getPairCorrelation(pairIdx);
    pair<unsigned int, unsigned int > layer;
    vector<double > minDistanceJ;
    minDistanceJ.resize(m_system.getNumOfType(tpair.second),1000.0);
    for (int atomTypeI=0; atomTypeI<m_system.getNumOfType(tpair.first); atomTypeI++)
      {
	layer.first=a_frame.getLayerOf(m_system.getIndexOfType(tpair.first, atomTypeI) );
	minDistance=1000.0;
	for (int atomTypeJ=0; atomTypeJ<m_system.getNumOfType(tpair.second); atomTypeJ++)
	  {
	    layer.second=a_frame.getLayerOf(m_system.getIndexOfType(tpair.second, atomTypeJ) );
	    pairDistance = a_frame.computeDistance(m_system.getIndexOfType(tpair.first,atomTypeI), m_system.getIndexOfType(tpair.second,atomTypeJ));
	    if (pairDistance < m_maxDist)
	      {
		binPairDistance(pairIdx, pairDistance);
		if (pairDistance < minDistanceJ[atomTypeJ])
		  {
		    minDistanceJ[atomTypeJ] = pairDistance;
		  }
		if (pairDistance < minDistance)
		  {
		    minDistance = pairDistance;
		  }
		incrementCounter(pairIdx);
	      }
	  }
	if(minDistance < m_maxDist)
          {
	    // bin distance as closest species J to species I (reference)
            binPairDistance(pairIdx, minDistance, 0, layer.first, layer.second);
          }
      }
    for (int atomTypeJ=0; atomTypeJ<m_system.getNumOfType(tpair.second); atomTypeJ++)
      {
	if(minDistanceJ[atomTypeJ] < m_maxDist)
          {
	    // bin distance as closest species I to species J (reference)
            binPairDistance(pairIdx, minDistanceJ[atomTypeJ], 1, layer.first, layer.second);
          }
      }
  }
}

void RDF::binPairDistance(unsigned int a_pair, double a_distance)
{
  int bin = floor(a_distance / m_binSize);
  m_rdf[a_pair][bin]++;
}

void RDF::binPairDistance(unsigned int a_pair, double a_distance, unsigned int a_firstLayer, unsigned int a_secondLayer)
{
  int bin = floor(a_distance / m_binSize);
  if (a_firstLayer == a_secondLayer )
    {
      m_rdfLayer[a_firstLayer][a_pair][bin]++;
    }
}

void RDF::binPairDistance(unsigned int a_pair, double a_distance, unsigned int a_whichClosest, unsigned int a_firstLayer, unsigned int a_secondLayer)
{
  int bin = floor(a_distance / m_binSize);
  if (a_firstLayer == a_secondLayer )
    {
      m_rdfLayerClosest[a_firstLayer][a_pair][a_whichClosest][bin]++;
    }
}

void RDF::incrementCounter(unsigned int a_pair)
{
  m_pairCounter[a_pair]++;
}

void RDF::normalize()
{
  for (int i=0; i<m_system.getNumPairs(); i++)
    {
      for (int j=0; j<m_numBins; j++)
	{
	  //m_rdf[i][j] /= m_pairCounter[i];
	  m_rdf[i][j] /= (4./3.)*M_PI*(pow(j+1,3)-pow(j,3))*pow(m_binSize,3);
	  for (int k=0; k<m_numLayers; k++)
	    {
	      m_rdfLayer[k][i][j] /= (4./3.)*M_PI*(pow(j+1,3)-pow(j,3))*pow(m_binSize,3);
	      m_rdfLayerClosest[k][i][0][j] /= (4./3.)*M_PI*(pow(j+1,3)-pow(j,3))*pow(m_binSize,3);
	      m_rdfLayerClosest[k][i][1][j] /= (4./3.)*M_PI*(pow(j+1,3)-pow(j,3))*pow(m_binSize,3);
	    }
	}
    }
}

void RDF::print()
{
  {
    for (int i=0; i<m_system.getNumPairs(); i++)
      {
	for (int j=0; j<m_numBins; j++)
	  {
	    cout << j*m_binSize << " " << m_rdf[i][j] << endl;
	  }
	cout << endl << endl;
      }
  }

}

const unsigned int RDF::getNumPairs() const
{
  return m_system.getNumPairs();
}
const unsigned int RDF::getNumBins() const
{
  return m_numBins;
}

const double RDF::getBinSize() const
{
  return m_binSize;
}

const double RDF::getRDFElement(int a_pair, int a_distance) const
{
  return m_rdf[a_pair][a_distance];
}

const double RDF::getRDFLayerElement(unsigned int a_layer, int a_pair, int a_distance) const
{
  return m_rdfLayer[a_layer][a_pair][a_distance];
}

const double RDF::getRDFLayerClosestElement(unsigned int a_layer, int a_pair, int a_whichClosest, int a_distance) const
{
  return m_rdfLayerClosest[a_layer][a_pair][a_whichClosest][a_distance];
}

const unsigned int RDF::getNumLayers() const
{
  return m_numLayers;
}

const char* RDFWrite(RDF* a_rdf, const char* a_filename)
{
  char full_filename [1024];

  if (strstr(a_filename, ".out") != NULL)
    {
      strcpy(full_filename, a_filename);
    }
  else
    {
      sprintf(full_filename, "%s.out", a_filename);
    }

  fp = fopen(full_filename, "w+");

  for (unsigned int iBin=0; iBin<a_rdf->getNumBins(); iBin++)
    {
      char str[128];
      sprintf(str, "%f", iBin*a_rdf->getBinSize());
      writeString(str);
      for (unsigned int jPair=0; jPair<a_rdf->getNumPairs(); jPair++)
	{
	  char str[128];
	  sprintf(str, " %f", a_rdf->getRDFElement(jPair,iBin));
	  writeString(str);
	}
      writeString("\n");
    }
  fclose(fp);

  return a_filename;
}

const char* RDFWriteLayers(RDF* a_rdf, const char* a_filename)
{
  unsigned int numLayers = a_rdf->getNumLayers();
  char** full_filenames = new char * [numLayers];

  if (strstr(a_filename, ".out") != NULL)
    {
      for (int iLayer=0; iLayer<numLayers; iLayer++)
	{
	  full_filenames[iLayer] = new char [1024];
	  strcpy(full_filenames[iLayer], a_filename);
	}
    }
  else
    {
      for (int iLayer=0; iLayer<numLayers; iLayer++)
	{
	  full_filenames[iLayer] = new char [1024];
	  sprintf(full_filenames[iLayer], "%s-%d.out", a_filename, iLayer);
	}
    }

  for (int iLayer=0; iLayer<numLayers; iLayer++)
    {
      fp = fopen(full_filenames[iLayer], "w+");
      writeString("#r ");
      for (unsigned int jPair=0; jPair<a_rdf->getNumPairs(); jPair++)
	{
	  char str[128];
	  sprintf(str, " %d_0 %d_1", jPair, jPair);
	  writeString(str);
	}
      writeString("\n");
      for (unsigned int iBin=0; iBin<a_rdf->getNumBins(); iBin++)
	{
	  char str[128];
	  sprintf(str, "%f", iBin*a_rdf->getBinSize());
	  writeString(str);
	  for (unsigned int jPair=0; jPair<a_rdf->getNumPairs(); jPair++)
	    {
	      char str[128];
	      sprintf(str, " %f %f", a_rdf->getRDFLayerClosestElement(iLayer,jPair, 0, iBin), a_rdf->getRDFLayerClosestElement(iLayer,jPair, 1, iBin));
	      writeString(str);
	    }
	  writeString("\n");
	}
      fclose(fp);
    }
  return a_filename;

}

