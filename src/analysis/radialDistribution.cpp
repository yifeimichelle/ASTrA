#include "frame.h"
#include "system.h"
#include "atom.h"
#include "radialDistribution.h"
#include <vector>
#include <iostream>
#include <assert.h>
#include "writer.h"

using namespace std;

RDF::RDF()
{
};

RDF::RDF(System& a_system)
{
  m_maxDist = 14.0;
  m_numBins = 500;
  m_system = a_system;
  m_numPairs = m_system.getNumPairs();
  m_numLayers = m_system.getNumLayers();
  m_rdf.resize(m_numBins);
  m_rdfLayer.resize(m_numLayers);
  m_rdfLayerClosest.resize(m_numLayers);

  m_pairCounter.resize(m_numPairs);

  for (int i=0; i<m_numBins; i++)
    {
      m_rdf[i].resize(m_numPairs);
    }
  for (int i=0; i<m_numLayers; i++)
    {
      m_rdfLayer[i].resize(m_numBins);
      m_rdfLayerClosest[i].resize(m_numBins);
      for (int j=0; j<m_numBins; j++)
	{
	  m_rdfLayer[i][j].resize(m_numPairs);
	  m_rdfLayerClosest[i][j].resize(m_numPairs);
	  for (int k=0; k<m_numPairs; k++)
	    {
	      m_rdfLayerClosest[i][j][k].resize(2);
	    }
	}
    }
  m_binSize = m_maxDist / m_numBins;
};

void RDF::sample(Frame& a_frame)
{
  double pairDistance;
  double minDistance;
  for (int pairIdx=0; pairIdx<m_numPairs; pairIdx++) {
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
		binPairDistance(pairDistance, pairIdx);
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
            binPairDistance(minDistance, pairIdx, 0, layer.first, layer.second);
          }
      }
    for (int atomTypeJ=0; atomTypeJ<m_system.getNumOfType(tpair.second); atomTypeJ++)
      {
	if(minDistanceJ[atomTypeJ] < m_maxDist)
          {
	    // bin distance as closest species I to species J (reference)
            binPairDistance(minDistanceJ[atomTypeJ], pairIdx, 1, layer.first, layer.second);
          }
      }
  }
}

void RDF::binPairDistance(double a_distance, unsigned int a_pair)
{
  int bin = floor(a_distance / m_binSize);
  m_rdf[bin][a_pair]++;
}

void RDF::binPairDistance(double a_distance, unsigned int a_pair, unsigned int a_firstLayer, unsigned int a_secondLayer)
{
  int bin = floor(a_distance / m_binSize);
  if (a_firstLayer == a_secondLayer )
    {
      m_rdfLayer[a_firstLayer][bin][a_pair]++;
    }
}

void RDF::binPairDistance(double a_distance, unsigned int a_pair, unsigned int a_whichClosest, unsigned int a_firstLayer, unsigned int a_secondLayer)
{
  int bin = floor(a_distance / m_binSize);
  if (a_firstLayer == a_secondLayer )
    {
      m_rdfLayerClosest[a_firstLayer][bin][a_pair][a_whichClosest]++;
    }
}

void RDF::incrementCounter(unsigned int a_pair)
{
  m_pairCounter[a_pair]++;
}

void RDF::normalize()
{
  for (int i=0; i<m_numBins; i++)
    {
      double normFactor = (4./3.)*M_PI*(pow(i+1,3)-pow(i,3))*pow(m_binSize,3);
      for (int j=0; j<m_system.getNumPairs(); j++)
	{
	  m_rdf[i][j] /= normFactor;
	  for (int k=0; k<m_numLayers; k++)
	    {
	      m_rdfLayer[k][i][j] /= normFactor;
	      m_rdfLayerClosest[k][i][j][0] /= normFactor;
	      m_rdfLayerClosest[k][i][j][1] /= normFactor;
	    }
	}
    }
}

void RDF::print()
{
  {
    for (int i=0; i<m_numBins; i++)
      {
	for (int j=0; j<m_system.getNumPairs(); j++)
	  {
	    cout << i*m_binSize << " " << m_rdf[i][j] << endl;
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

const double RDF::getRDFElement(int a_distance, int a_pair) const
{
  return m_rdf[a_distance][a_pair];
}

const double RDF::getRDFLayerElement(unsigned int a_layer, int a_distance, int a_pair) const
{
  return m_rdfLayer[a_layer][a_distance][a_pair];
}

const double RDF::getRDFLayerClosestElement(unsigned int a_layer, int a_distance, int a_pair, int a_whichClosest) const
{
  return m_rdfLayerClosest[a_layer][a_distance][a_pair][a_whichClosest];
}

const unsigned int RDF::getNumLayers() const
{
  return m_numLayers;
}

double* RDF::getRDFAddress(int i)
{
  return &(m_rdf[i][0]);
}

void RDF::setRDFLayerClosestValue(int a_layer, int a_bin, int a_pair, int a_closest, double a_setVal)
{
    m_rdfLayerClosest[a_layer][a_bin][a_pair][a_closest] = a_setVal;
}


double*** RDF::getRDFAddressLayers(int i)
{
  return (double***)&(m_rdfLayerClosest[i][0]); // how to typecast to pointer to pointer?
}

double** RDF::getRDFAddressLayers(int i, int j)
{
  return (double**)&(m_rdfLayerClosest[i][j][0]);
}

double* RDF::getRDFAddressLayers(int i, int j, int k)
{
  return &(m_rdfLayerClosest[i][j][k][0]);
}


const char* RDFWrite(RDF* a_rdf, const char* a_filename)
{
  double binSize = a_rdf->getBinSize();
  int numBins = a_rdf->getNumBins();
  int varDim = a_rdf->getNumPairs();
  const char * const headernames[] = { "nodeData" };
  double* data[500];
  for (int i=0; i<numBins; i++)
    {
      data[i] = a_rdf->getRDFAddress(i);
    }
  write_binned_data(a_filename, numBins, binSize, varDim, headernames, data);
  return a_filename;
}

const char* RDFWriteLayers(RDF* a_rdf, const char* a_filename)
{
  double binSize = a_rdf->getBinSize();
  int numBins = a_rdf->getNumBins();
  int varDim = a_rdf->getNumPairs();
  int numLayers = a_rdf->getNumLayers();
  const char * const headernames[] = { "nodeData" };
  double**** data = new double***[numLayers];
  for (int i=0; i<numLayers; i++)
    {
        data[i] = new double**[numBins];
      for (int j=0; j<numBins; j++)
      	{
            data[i][j] = new double*[varDim];
          for (int k=0; k<varDim; k++)
          {
            data[i][j][k] = a_rdf->getRDFAddressLayers(i,j,k);
          }
      	}
    }
  write_binned_layered_data(a_filename, numBins, binSize, varDim, numLayers, 2, headernames, data);

  // delete array of pointers
  for (int i=0; i<numLayers; i++)
    {
      for (int j=0; j<numBins; j++)
      	{
            delete data[i][j];
      	}
    }
  for (int i=0; i<numLayers; i++)
    {
        delete data[i];
    }
  delete data;

  return a_filename;
}

