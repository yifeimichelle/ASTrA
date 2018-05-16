#include "frame.h"
#include "system.h"
#include "atom.h"
#include "radialDistribution.h"
#include <vector>
#include <iostream>
#include <assert.h>
#include "writer.h"
#include <math.h>

using namespace std;

RDF::RDF()
{
};

RDF::RDF(System& a_system)
{
  m_maxDist = 14.0;
  m_numBins = 500;
  m_binSize = m_maxDist / m_numBins;
  m_system = a_system;
  m_numPairs = m_system.getNumPairs();
  m_numMolecPairs = m_system.getNumMolecPairs();
  m_numLayers = m_system.getNumLayers();
  m_rdf.resize(m_numBins);
  m_rdfLayer.resize(m_numLayers);
  m_rdfLayerClosest.resize(m_numLayers);
  m_rdfMolec.resize(m_numBins);
  m_rdfMolecLayer.resize(m_numLayers);
  m_currentRDFMolecLayer.resize(m_numLayers);
  m_rdfMolecLayerClosest.resize(m_numLayers);

  // Constants for calculating DoC
  m_Rj = 0.715; // half of carbon-carbon distance in SP2 structure
  m_phi = 0.6046; // coverage of surface covered by hexagonally tiled structure
  m_Rcut = 6.3;
  m_solidAngleFactor.resize(m_numBins);
  m_numBinsDoC = 100;
  m_binSizeDoC = 1.0/m_numBinsDoC;
  double RjSq = m_Rj * m_Rj;
  for (int i=0; i<m_numBins; i++)
    {
      double rBin = i*m_binSize;
      double factor;
      factor = rBin*rBin + RjSq;
      factor = sqrt (factor);
      factor = rBin / factor;
      m_solidAngleFactor[i] = 1 - factor;
    }
  m_DoC.resize(m_numBins);
  for (int i=0; i<m_numBins; i++)
    {
      m_DoC[i].resize(m_numMolecPairs);
    }
  m_DoCHist.resize(m_numBinsDoC);
  for (int i=0; i<m_numBinsDoC; i++)
    {
      m_DoCHist[i].resize(m_numMolecPairs);
    }
  // end constants for calculating DoC

  m_pairCounter.resize(m_numPairs);
  m_pairMolecCounter.resize(m_numMolecPairs);

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

  for (int i=0; i<m_numBins; i++)
    {
      m_rdfMolec[i].resize(m_numMolecPairs);
    }
  for (int i=0; i<m_numLayers; i++)
    {
      m_rdfMolecLayer[i].resize(m_numBins);
      m_currentRDFMolecLayer[i].resize(m_numBins);
      m_rdfMolecLayerClosest[i].resize(m_numBins);
      for (int j=0; j<m_numBins; j++)
	{
	  m_rdfMolecLayer[i][j].resize(m_numMolecPairs);
	  m_currentRDFMolecLayer[i][j].resize(m_numMolecPairs);
	  m_rdfMolecLayerClosest[i][j].resize(m_numMolecPairs);
	  for (int k=0; k<m_numMolecPairs; k++)
	    {
	      m_rdfMolecLayerClosest[i][j][k].resize(2);
	    }
	}
    }

};

void RDF::sampleZP(Frame& a_frame)
{
  // FIXME : Michelle
  // Do I need this routine?
}

void RDF::sample(const Frame& a_frame)
{
  sampleMolecules(a_frame);
  sampleAtoms(a_frame);
}

void RDF::sampleMolecules(const Frame& a_frame)
{
  double pairDistance;
  double minDistance;
  // For each pair
  //   For each layer
  //     Iterate through atoms (or ions) in layer
  //       Compute distance and bin

  // For each layer
  for (int layIdx = 0; layIdx < m_numLayers; layIdx++)
    {
      // Get list of molecules in this layer
      vector<int>* molecsInLayer = a_frame.getMoleculesInLayer(layIdx);
      // For each pair
      for (int pairIdx = 0; pairIdx < m_numMolecPairs; pairIdx++)
	{
	  pair<unsigned int, unsigned int > molecPair = m_system.getMolecPairCorrelation(pairIdx);
	  int pairFirst = molecPair.first;
	  int pairSecond = molecPair.second;
	  vector<double > minDistanceB;
	  int pairSecondSize = molecsInLayer[pairSecond].size();
	  minDistanceB.resize(pairSecondSize,1000.0);

	  // Find out whether to compute DoC
	  int layer = -1;
	  unsigned int computeDoC = 0;
	  if ( m_system.isCathode( pairSecond ) or  m_system.isCathode( pairFirst ))
	    {
	      if ( m_system.isAnodeLower() )
		{
		  layer = 2;
		}
	      else
		{
		  layer = 0;
		}
	    }
	  else if ( m_system.isAnode( pairSecond ) or m_system.isAnode( pairFirst ) )
	    {
	      if ( m_system.isAnodeLower() )
		{
		  layer = 0;
		}
	      else
		{
		  layer = 2;
		}
	    }
	  if ( layIdx == layer)
	    {
	      computeDoC = 1;
	    }
		
	  // For each molecule of first type in pair
	  for (vector<int>::iterator itA = molecsInLayer[pairFirst].begin(); itA != molecsInLayer[pairFirst].end(); ++itA)
	    {
	      minDistance = 1000.0;
	      int secondIndex = 0;
	      double DoC = 0;
	      // For each molecule of second type in pair
	      for (vector<int>::iterator itB = molecsInLayer[pairSecond].begin(); itB != molecsInLayer[pairSecond].end(); ++itB)
		{
		  // Compute distance
		  float distance = a_frame.computeMolecDistance(*itA,*itB);
		  // Check whether distance is closest for atom type A
		  if (distance < minDistance)
		    {
		      minDistance = distance;
		    }
  		  // Check whether distance is closest for atom type B
		  if (distance < minDistanceB[secondIndex])
		    {
		      minDistanceB[secondIndex] = distance;
		    }
		  // Bin it
		  if (distance < m_maxDist)
		    {
		      binMolecPairDistanceLayer(distance, pairIdx, layIdx);
		      if (computeDoC)
			{
			  if (distance < m_Rcut)
			    {
			      DoC += computeSolidAngleFactor(distance);
			    }
			}

		    }
		  secondIndex++;
		}
	      DoC /= (2.0*m_phi);
	      binDoC(DoC, pairIdx);
	      if(minDistance < m_maxDist)
	      	{
	      	  // bin distance as closest species B to species A (reference)
	      	  binMolecPairDistanceClosestLayer(minDistance, pairIdx, 0, layIdx);
	      	}
	    }
	  // For each molecule of second type in pair
	  for (int secondIndex=0; secondIndex<pairSecondSize; secondIndex++)
	    {
	      if(minDistanceB[secondIndex] < m_maxDist)
	  	{
	  	  // bin distance as closest species I to species J (reference)
	  	  binMolecPairDistanceClosestLayer(minDistanceB[secondIndex], pairIdx, 1, layIdx);
	  	}
	    }
	}
    }

}

void RDF::sampleAtoms(const Frame& a_frame)
{
  double pairDistance;
  double minDistance;
  // For each pair
  //   For each layer
  //     Iterate through atoms (or ions) in layer
  //       Compute distance and bin

  // For each layer
  for (int layIdx = 0; layIdx < m_numLayers; layIdx++)
    {
      // Get list of atoms in this layer
      vector<int>* atomsInLayer = a_frame.getAtomsInLayer(layIdx);
      // For each pair
      for (int pairIdx = 0; pairIdx < m_numPairs; pairIdx++)
	{
	  // Get atom types of pair
	  pair<unsigned int, unsigned int > atomPair = m_system.getPairCorrelation(pairIdx);
	  int pairFirst = atomPair.first;
	  int pairSecond = atomPair.second;
	  int pairSecondSize = atomsInLayer[pairSecond].size();
	  vector<double > minDistanceB;
	  minDistanceB.resize(pairSecondSize,1000.0);
	  // For each atom of first type in pair
	  for (vector<int>::iterator itA = atomsInLayer[pairFirst].begin(); itA != atomsInLayer[pairFirst].end(); ++itA)
	    {
	      minDistance = 1000.0;
	      int secondIndex = 0;
	      // For each atom of second type in pair
	      for (vector<int>::iterator itB = atomsInLayer[pairSecond].begin(); itB != atomsInLayer[pairSecond].end(); ++itB)
		{
		  // Compute distance
		  float distance = a_frame.computeDistance(*itA,*itB);
		  // Check whether distance is closest for atom type A
		  if (distance < minDistance)
		    {
		      minDistance = distance;
		    }
		  // Check whether distance is closest for atom type B
		  if (distance < minDistanceB[secondIndex])
		    {
		      minDistanceB[secondIndex] = distance;
		    }
		  // Bin it
		  if (distance < m_maxDist)
		    {
		      binPairDistanceLayer(distance, pairIdx, layIdx);
		    }
		  secondIndex++;
		}
	      if(minDistance < m_maxDist)
	      	{
	      	  // bin distance as closest species B to species A (reference)
	      	  binPairDistanceClosestLayer(minDistance, pairIdx, 0, layIdx);
	      	}
	    }
	  // For each atom of second type in pair
	  for (int secondIndex=0; secondIndex<pairSecondSize; secondIndex++)
	    {
	      if(minDistanceB[secondIndex] < m_maxDist)
	  	{
	  	  // bin distance as closest species I to species J (reference)
	  	  binPairDistanceClosestLayer(minDistanceB[secondIndex], pairIdx, 1, layIdx);
	  	}
	    }
	}
    }

}

void RDF::binPairDistanceLayer(double a_distance, unsigned int a_pair, unsigned int a_layer)
{
  int bin = floor(a_distance / m_binSize);
  m_rdf[bin][a_pair]++;
  m_rdfLayer[a_layer][bin][a_pair]++;
}

void RDF::binPairDistanceClosestLayer(double a_distance, unsigned int a_pair, unsigned int a_whichClosest, unsigned int a_layer)
{
  int bin = floor(a_distance / m_binSize);
  m_rdfLayerClosest[a_layer][bin][a_pair][a_whichClosest]++;
}

void RDF::binMolecPairDistanceLayer(double a_distance, unsigned int a_pair, unsigned int a_layer)
{
  int bin = floor(a_distance / m_binSize);
  m_rdfMolec[bin][a_pair]++;
  m_rdfMolecLayer[a_layer][bin][a_pair]++;
  m_currentRDFMolecLayer[a_layer][bin][a_pair]++;
}

void RDF::clearFrame()
{
    for (int i=0; i<m_numLayers; i++)
    {
      for (int j=0; j<m_numBins; j++)
	{
	  fill(m_currentRDFMolecLayer[i][j].begin(), m_currentRDFMolecLayer[i][j].end(), 0);
	}
    }
}

void RDF::binMolecPairDistanceClosestLayer(double a_distance, unsigned int a_pair, unsigned int a_whichClosest, unsigned int a_layer)
{
  int bin = floor(a_distance / m_binSize);
  m_rdfMolecLayerClosest[a_layer][bin][a_pair][a_whichClosest]++;
}

void RDF::incrementCounter(unsigned int a_pair)
{
  m_pairCounter[a_pair]++;
}

void RDF::incrementMolecCounter(unsigned int a_pair)
{
  m_pairMolecCounter[a_pair]++;
}

double RDF::computeSolidAngleFactor(double a_distance)
{
  double factor;
  double retVal;
  double RjSq = m_Rj * m_Rj;
  factor = a_distance*a_distance + RjSq;
  factor = sqrt (factor);
  factor = a_distance / factor;
  retVal = 1.0 - factor;
  return retVal;
}

double RDF::binDoC(double a_doc, unsigned int a_pair)
{
  int bin = floor(a_doc / m_binSizeDoC);
  assert(a_doc<1.0);
  m_DoCHist[bin][a_pair]++;
  m_countIonsDoC[a_pair]++;
}


// Computes degree of confinement after number of ions in rdf have been counted
void RDF::computeDegreeOfConfinement(const Frame& a_frame)
{
  // For molec pair
  for (int pairIdx = 0; pairIdx < m_numMolecPairs; pairIdx++)
    {
      //cout << "computing DoC" << endl;
      //cout << pairIdx << endl;
      pair<unsigned int, unsigned int > atomPair = m_system.getMolecPairCorrelation(pairIdx);
      unsigned int pairFirst = atomPair.first;
      unsigned int pairSecond = atomPair.second;
      int layer = -1;
      double sum = 0.0;
      // Layer = bottom if 2nd type is anode, top if cathode, skip if neither
      if ( m_system.isCathode( pairSecond ) or  m_system.isCathode( pairFirst ))
      	{
	  if ( m_system.isAnodeLower() )
	    {
	      layer = 2;
	    }
	  else
	    {
	      layer = 0;
	    }
      	}
      else if ( m_system.isAnode( pairSecond ) or m_system.isAnode( pairFirst ) )
      	{
	  if ( m_system.isAnodeLower() )
	    {
	      layer = 0;
	    }
	  else
	    {
	      layer = 2;
	    }
      	}
      // If molecule is anode or cathode, layer will be set to something nonzero. Otherwise don't compute degree of confinement.
      if (layer > -1)
	{
	  unsigned int numIons = a_frame.getCurrentNumMolecsInLayer(layer, pairFirst);
	  for (int binIdx = 0; binIdx < m_numBins; binIdx ++)
	    {
	      // increment interior sum
	      sum += m_solidAngleFactor[binIdx]*m_currentRDFMolecLayer[layer][binIdx][pairIdx];
	      // save interior sum
	      m_DoC[binIdx][pairIdx] += sum/(numIons);
	      //cout << "adding to DoC ..." << m_currentRDFMolecLayer[layer][binIdx][pairIdx] << " "< < numIons << endl;
	    }
	}
    }
}

/// Get number of bins over which the DoC is computed.
const unsigned int RDF::getNumBinsDoC() const
{
  return m_numBinsDoC;
}
/// Get size of bins of DoC (in length units).
const double RDF::getBinSizeDoC() const
{
  return m_binSizeDoC;
}

void RDF::normalize()
{
  int numFrames = m_system.getNumFrames();
  for (int i=0; i<m_numBins; i++)
    {
      double normFactor = (4./3.)*M_PI*(pow(i+1,3)-pow(i,3))*pow(m_binSize,3);
      normFactor = normFactor * numFrames;
      double DoCnormFactor = numFrames * m_phi*2.0 * m_binSizeDoC;
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
      for (int j=0; j<m_system.getNumMolecPairs(); j++)
	{
	  m_rdfMolec[i][j] /= normFactor;
	  m_DoC[i][j] /= DoCnormFactor;
	  for (int k=0; k<m_numLayers; k++)
	    {
	      m_rdfMolecLayer[k][i][j] /= normFactor;
	      m_rdfMolecLayerClosest[k][i][j][0] /= normFactor;
	      m_rdfMolecLayerClosest[k][i][j][1] /= normFactor;
	    }
	}
    }
  for (int i=0; i<m_numBinsDoC; i++)
    {
      for (int j=0; j<m_system.getNumMolecPairs(); j++)
	{
	  m_DoCHist[i][j] /= (numFrames * m_countIonsDoC[j]);
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

const unsigned int RDF::getNumMolecPairs() const
{
  return m_system.getNumMolecPairs();
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

double* RDF::getRDFAddressLayers(int i, int j)
{
  return &(m_rdfLayer[i][j][0]);
}

double* RDF::getRDFAddressLayersClosest(int i, int j, int k)
{
  return &(m_rdfLayerClosest[i][j][k][0]);
}

double* RDF::getMolecRDFAddress(int i)
{
  return &(m_rdfMolec[i][0]);
}

double* RDF::getMolecRDFAddressLayers(int i, int j)
{
  return &(m_rdfMolecLayer[i][j][0]);
}

double* RDF::getMolecRDFAddressLayersClosest(int i, int j, int k)
{
  return &(m_rdfMolecLayerClosest[i][j][k][0]);
}

double* RDF::getDoCAddress(int i)
{
  return &(m_DoC[i][0]);
}

double* RDF::getDoCHistAddress(int i)
{
  return &(m_DoCHist[i][0]);
}

const char* RDFWrite(RDF* a_rdf, const char* a_filename)
{
  double binSize = a_rdf->getBinSize();
  int numBins = a_rdf->getNumBins();
  int varDim = a_rdf->getNumPairs();
  const char * const headernames[] = { "z[A]",  "0",  "1",  "2",  "3",  "4",  "5",  "6",  "7",  "8",  "9",  "10",  "11",  "12",  "13", "14", "15", "16", "17", "18", "19", "20" };
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
  double*** data = new double**[numLayers];
  for (int i=0; i<numLayers; i++)
    {
      data[i] = new double*[numBins];
      for (int j=0; j<numBins; j++)
      	{
	  data[i][j] =a_rdf->getRDFAddressLayers(i,j);
	}
    }
  write_binned_layered_data(a_filename, numBins, binSize, varDim, numLayers, headernames, data);

  // delete array of pointers
  for (int i=0; i<numLayers; i++)
    {
        delete data[i];
    }
  delete data;

  return a_filename;
}

const char* RDFWriteLayersClosest(RDF* a_rdf, const char* a_filename)
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
	      data[i][j][k] = a_rdf->getRDFAddressLayersClosest(i,j,k);
	    }
      	}
    }
  write_binned_layered_multival_data(a_filename, numBins, binSize, varDim, numLayers, 2, headernames, data);

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


const char* RDFMolecWrite(RDF* a_rdf, const char* a_filename)
{
  double binSize = a_rdf->getBinSize();
  int numBins = a_rdf->getNumBins();
  int varDim = a_rdf->getNumMolecPairs();
  const char * const headernames[] = { "z[A]",  "0",  "1",  "2",  "3",  "4",  "5",  "6",  "7",  "8",  "9",  "10",  "11",  "12",  "13", "14", "15", "16", "17", "18", "19", "20" };
  double* data[500];
  for (int i=0; i<numBins; i++)
    {
      data[i] = a_rdf->getMolecRDFAddress(i);
    }
  write_binned_data(a_filename, numBins, binSize, varDim, headernames, data);
  return a_filename;
}

const char* RDFMolecWriteLayers(RDF* a_rdf, const char* a_filename)
{
  double binSize = a_rdf->getBinSize();
  int numBins = a_rdf->getNumBins();
  int varDim = a_rdf->getNumMolecPairs();
  int numLayers = a_rdf->getNumLayers();
  const char * const headernames[] = { "nodeData" };
  double*** data = new double**[numLayers];
  for (int i=0; i<numLayers; i++)
    {
      data[i] = new double*[numBins];
      for (int j=0; j<numBins; j++)
      	{
	  data[i][j] =a_rdf->getMolecRDFAddressLayers(i,j);
	}
    }
  write_binned_layered_data(a_filename, numBins, binSize, varDim, numLayers, headernames, data);

  // delete array of pointers
  for (int i=0; i<numLayers; i++)
    {
        delete data[i];
    }
  delete data;

  return a_filename;
}


const char* RDFMolecWriteLayersClosest(RDF* a_rdf, const char* a_filename)
{
  double binSize = a_rdf->getBinSize();
  int numBins = a_rdf->getNumBins();
  int varDim = a_rdf->getNumMolecPairs();
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
	      data[i][j][k] = a_rdf->getMolecRDFAddressLayersClosest(i,j,k);
	    }
      	}
    }
  write_binned_layered_multival_data(a_filename, numBins, binSize, varDim, numLayers, 2, headernames, data);

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

const char* DoCWrite(RDF* a_rdf, const char* a_filename)
{
  double binSize = a_rdf->getBinSize();
  int numBins = a_rdf->getNumBins();
  int varDim = a_rdf->getNumMolecPairs();
  const char * const headernames[] = { "z[A]",  "0",  "1",  "2",  "3",  "4",  "5",  "6",  "7",  "8",  "9",  "10",  "11",  "12",  "13", "14", "15", "16", "17", "18", "19", "20" };
  double* data[500];
  for (int i=0; i<numBins; i++)
    {
      data[i] = a_rdf->getDoCAddress(i);
    }
  write_binned_data(a_filename, numBins, binSize, varDim, headernames, data);
  return a_filename;
}

const char* DoCHistWrite(RDF* a_rdf, const char* a_filename)
{
  double binSize = a_rdf->getBinSizeDoC();
  int numBins = a_rdf->getNumBinsDoC();
  int varDim = a_rdf->getNumMolecPairs();
  const char * const headernames[] = { "z[A]",  "0",  "1",  "2",  "3",  "4",  "5",  "6",  "7",  "8",  "9",  "10",  "11",  "12",  "13", "14", "15", "16", "17", "18", "19", "20" };
  double* data[500];
  for (int i=0; i<numBins; i++)
    {
      data[i] = a_rdf->getDoCHistAddress(i);
    }
  write_binned_data(a_filename, numBins, binSize, varDim, headernames, data);
  return a_filename;
}
