#include "frame.h"
#include "system.h"
#include "atomCounter.h"
#include "atom.h"
#include "radialDistribution.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <assert.h>
#include "writer.h"
#include <math.h>

using namespace std;

RDF::RDF()
{
};

RDF::~RDF()
{
  for (int i=0; i<m_system.getNumDoCThresholds(); i++)
  {
    m_DoCIndicesFile[i].close();
  }
}

RDF::RDF(System& a_system, AtomCounter& a_ac)
{
  m_DoCIndicesFile.resize(a_system.getNumDoCThresholds());
  for (int i=0; i<a_system.getNumDoCThresholds(); i++)
  {
    char filename[1024];
    char full_filename[1024];
    sprintf(filename, "%s-%d", "DoCIndices", i);
    sprintf(full_filename, "%s.out", filename);
    m_DoCIndicesFile[i].open(full_filename, ios::out | ios::trunc);
  }
  m_maxDist = 14.0;
  m_numBins = 500;
  m_binSize = m_maxDist / m_numBins;
  m_system = a_system;
  m_ac = a_ac;
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

  // For calculating coordination number
  m_numBinsCoordNum = MAX_COORD_NUM;
  m_binSizeCoordNum = MAX_COORD_NUM * 1.0 / m_numBinsCoordNum;
  m_coordNum.resize(m_numLayers);
  m_coordNumHist.resize(m_numLayers);
  m_countCoordNum.resize(m_numLayers);
  for (int i=0; i<m_numLayers; i++)
  {
    m_coordNum[i].resize(m_numMolecPairs);
    m_countCoordNum[i].resize(m_numMolecPairs);
    m_coordNumHist[i].resize(m_numBinsCoordNum);
    for (int j=0; j<m_numBinsCoordNum; j++)
    {
      m_coordNumHist[i][j].resize(m_numMolecPairs);
    }
  }

  // Constants for calculating DoC
  m_Rj = 0.715; // half of carbon-carbon distance in SP2 structure
  m_phi = 0.6046; // coverage of surface covered by hexagonally tiled structure
  m_ionCharge = 0.78;
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
    m_DoCHist[i].resize(3*m_numMolecPairs);
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

const System& RDF::getSystem() const
{
  return m_system;
}

const AtomCounter& RDF::getAtomCounter() const
{
  return m_ac;
}

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
      double Rcut = m_system.getMolecPairCutoff(pairIdx);
      int pairFirst = molecPair.first;
      int pairSecond = molecPair.second;
      vector<double > minDistanceB;
      int pairSecondSize = molecsInLayer[pairSecond].size();
      minDistanceB.resize(pairSecondSize,1000.0);

      // Find out whether to compute DoC
      int layer = -1;
      unsigned int computeDoC = 0;
      unsigned int isCounterCharge = 0;
      // If any of the atoms in the pair are electrode atoms, then rdf/doc only exists in that electrode's layer
      if ( m_system.isCathode( pairSecond ) or  m_system.isCathode( pairFirst ))
      {
        if ( m_system.isCathodeLower() )
        {
          layer = 0;
        }
        else
        {
          layer = 2;
        }
        if ( m_system.isCation(pairFirst ) )
        {
          isCounterCharge = 1;
        }
      }
      else if ( m_system.isAnode( pairSecond ) or m_system.isAnode( pairFirst ) )
      {
        if ( m_system.isCathodeLower() )
        {
          layer = 2;
        }
        else
        {
          layer = 0;
        }
        if ( m_system.isAnion(pairFirst ) )
        {
          isCounterCharge = 1;
        }
      }
      if ( layIdx == layer)
      {
        computeDoC = 1;
      }
      // For each molecule of first type in pair
      //cout << layIdx << " " << pairIdx << " " << molecsInLayer[pairFirst].size() << " " << molecsInLayer[pairSecond].size() << endl;
      for (vector<int>::iterator itA = molecsInLayer[pairFirst].begin(); itA != molecsInLayer[pairFirst].end(); ++itA)
      {
        minDistance = 1000.0;
        int secondIndex = 0;
        double DoC = 0.0;
        double sumElecCharge = 0.0;
        int numCoordCarbons = 0;
        double coordNum = 0.0;
        // For each molecule of second type in pair
        for (vector<int>::iterator itB = molecsInLayer[pairSecond].begin(); itB != molecsInLayer[pairSecond].end(); ++itB)
        {
          if (*itA != *itB) // exclude self-self pair
          {
            // Compute distance
            double distance = a_frame.computeMolecDistance(*itA,*itB);
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
              if (distance < Rcut)
              {
                if (computeDoC)
                {
                  DoC += computeSolidAngleFactor(distance);
                  sumElecCharge += a_frame.getAtomOfMolec(*itB).getCharge();
                  numCoordCarbons += 1;
                }
                else
                {
                  coordNum++;
                  // if (coordNum > 20)
                  // {
                  //   cout << pairIdx << " " <<  m_system.getFirstAtomOfMolec(*itA);
                  // 	 cout << " " << m_system.getFirstAtomOfMolec(*itB) << " ";
                  // 	 cout << coordNum << " " << Rcut << " " << distance << endl;
                  // }
                }
              }
            }
          }
          // increment index for minDistanceB
          secondIndex++;
        }
        if(coordNum > MAX_COORD_NUM)
        {
          cout << coordNum << endl;
        }
        binCoordNum(coordNum, pairIdx, layIdx);
        if (computeDoC)
        {
          DoC /= (2.0*m_phi);
          binDoC(a_frame, m_system.getFirstAtomOfMolec(*itA), DoC, sumElecCharge, isCounterCharge, pairIdx, numCoordCarbons);
        }
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

double RDF::binDoC(const Frame& a_frame, int a_atomID, double a_doc, double a_elecCharge, unsigned int a_isCounterCharge, unsigned int a_pair, int a_numCoordCarbons)
{
  int bin = floor(a_doc / m_binSizeDoC);
  if(a_doc<1.0)
  {
    m_DoCHist[bin][3*a_pair]++;
    m_DoCHist[bin][3*a_pair+1] += a_elecCharge / m_ionCharge;
    m_DoCHist[bin][3*a_pair+2] += a_numCoordCarbons;
    m_countIonsDoC[a_pair]++;
    for (int i=0; i<m_system.getNumDoCThresholds(); i++)
    {
      if (a_doc > m_system.getDoCThresholdLo(i))
      {
        if (a_doc < m_system.getDoCThresholdHi(i))
        {
          const Atom& atom = a_frame.getAtom(a_atomID);
          const array<double, DIM>& pos = atom.getPosition();

          // save timestep and atom indices
          m_DoCIndicesFile[i] << " " << atom.getName() << " " << pos[0] << " " << pos[1] << " " << pos[2];
          m_DoCIndicesFile[i] << " " << a_doc << " " << a_elecCharge << " " << a_numCoordCarbons;
          m_DoCIndicesFile[i] << " " << a_atomID << " " << a_frame.getTimestep() << endl;
          //writeDoCIndicesToFile();
        }
      }
    }
  }
}

double RDF::binCoordNum(double a_coordNum, unsigned int a_pair, unsigned int a_layer)
{
  int bin = floor(a_coordNum / m_binSizeCoordNum);
  // if(a_coordNum < MAX_COORD_NUM)
  //   {
  assert(a_coordNum < MAX_COORD_NUM);
  m_coordNumHist[a_layer][bin][a_pair]++;
  m_coordNum[a_layer][a_pair] += a_coordNum;
  m_countCoordNum[a_layer][a_pair]++;
  //   }
  // else if (a_coordNum > MAX_COORD_NUM)
  //   {
  //      if (m_errorPrinted == 0)
  // 	{
  // 	  cout << "WARNING: Some coordination numbers are greater than " << MAX_COORD_NUM << "!" << endl;
  // 	  m_errorPrinted = 1;
  // 	}
  //   }
}


// Computes degree of confinement after number of ions in rdf have been counted
void RDF::computeDegreeOfConfinement(const Frame& a_frame)
{
  //cout << "starting to compute degree of confinement." << endl;
  // For molec pair
  for (int pairIdx = 0; pairIdx < m_numMolecPairs; pairIdx++)
  {
    pair<unsigned int, unsigned int > molecPair = m_system.getMolecPairCorrelation(pairIdx);
    unsigned int pairFirst = molecPair.first;
    unsigned int pairSecond = molecPair.second;
    int layer = -1;
    double sum = 0.0;
    // Layer = bottom if 2nd type is anode, top if cathode, skip if neither
    if ( m_system.isCathode( pairSecond ) or  m_system.isCathode( pairFirst ))
    {
      if ( m_system.isCathodeLower() )
      {
        layer = 0;
      }
      else
      {
        layer = 2;
      }
    }
    else if ( m_system.isAnode( pairSecond ) or m_system.isAnode( pairFirst ) )
    {
      if ( m_system.isCathodeLower() )
      {
        layer = 2;
      }
      else
      {
        layer = 0;
      }
    }
    // If molecule is anode or cathode, layer will be set to something nonzero. Otherwise don't compute degree of confinement.
    if (layer > -1)
    {
      //cout << "computing degree of confinement for pair " << pairIdx <<" in layer " << layer << "." << endl;
      unsigned int numIons = a_frame.getCurrentNumMolecsInLayer(layer, pairFirst);
      for (int binIdx = 0; binIdx < m_numBins; binIdx ++)
      {
        //cout << "going through bin " << binIdx << " of " << m_numBins << " for pair " << pairIdx <<" in layer " << layer << "." << endl;
        // increment interior sum
        sum += m_solidAngleFactor[binIdx]*m_currentRDFMolecLayer[layer][binIdx][pairIdx];
        // save interior sum
        m_DoC[binIdx][pairIdx] += sum/(numIons);
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


const unsigned int RDF::getNumBinsCoordNum() const
{
  return m_numBinsCoordNum;
}

const double RDF::getBinSizeCoordNum() const
{
  return m_binSizeCoordNum;
}


void RDF::normalize(AtomCounter* a_ac)
{
  // Store layer upper bounds in volLayer
  double* volLayer;
  volLayer = new double[m_numLayers];
  m_system.getLayerUpperBounds(m_numLayers,volLayer);
  // Subtract to get height of layer
  for (int i=m_numLayers-1; i>0; i--)
  {
    volLayer[i] = volLayer[i] - volLayer[i-1];
  }
  // Get cross sectional area of unit cell
  double boxArea = m_system.getBoxDim(0)*m_system.getBoxDim(1);
  // Create array to hold number of ions in layer
  double** avgIonsInLayer;
  avgIonsInLayer = new double* [m_numLayers];
  // Multiply layer height by area to get layer volume
  // Assign pointers to address of ions in layer from AtomCounter
  for (int i=0; i<m_numLayers; i++)
  {
    volLayer[i] *= boxArea;
    avgIonsInLayer[i] = a_ac->getACIonsLayersAddress(i);
  }

  int numFrames = m_system.getNumFrames();
  // For each bin
  for (int i=0; i<m_numBins; i++)
  {
    // contribution of volume to normalization factor
    // (4/3)*pi*(R^3-r^3)
    double normFactor = (4./3.)*M_PI*(pow(i+1,3)-pow(i,3))*pow(m_binSize,3);
    normFactor = normFactor * numFrames;
    // normalization for degree of confinement
    double DoCnormFactor = numFrames * m_phi * 2.0;
    // For each atom-atom pair
    for (int j=0; j<m_system.getNumPairs(); j++)
    {
      // Normalize overall rdf
      m_rdf[i][j] /= normFactor;
      for (int k=0; k<m_numLayers; k++)
      {
        // Normalize per-layer rdf
        m_rdfLayer[k][i][j] /= normFactor;
        m_rdfLayerClosest[k][i][j][0] /= normFactor;
        m_rdfLayerClosest[k][i][j][1] /= normFactor;
      }
    }
    // For each molecule-molecule pair
    for (int j=0; j<m_numMolecPairs; j++)
    {
      // Get molecule IDs for pair
      pair<unsigned int, unsigned int > molecPair = m_system.getMolecPairCorrelation(j);
      int pairFirst = molecPair.first;
      int pairSecond = molecPair.second;

      // Normalize overall rdf and DoC
      m_rdfMolec[i][j] /= normFactor;
      m_DoC[i][j] /= DoCnormFactor;

      // Get layer- and ion-specific normalization factor for rdfMolecLayer
      for (int k=0; k<m_numLayers; k++)
      {
        double densNormFactor;
        if (pairFirst < 3 && pairSecond < 3) {
          // Only do this for cation, anion, and solvent (molecules 0, 1, 2)
          densNormFactor = avgIonsInLayer[k][pairFirst]*avgIonsInLayer[k][pairSecond] / volLayer[k];
          //cout << k << " " << j << " " << pairFirst << " " << pairSecond << " ";
          //cout << avgIonsInLayer[k][pairFirst] << " " << avgIonsInLayer[k][pairSecond] << " " << volLayer[k] << " " << densNormFactor << endl;
        }
        else
        {
          // Otherwise don't multiply by factor
          densNormFactor=1;
        }
        // Normalize per-layer rdf
        m_rdfMolecLayer[k][i][j] /= normFactor * densNormFactor;
        m_rdfMolecLayerClosest[k][i][j][0] /= normFactor;
        m_rdfMolecLayerClosest[k][i][j][1] /= normFactor;
      }
    }
  }
  for (int i=0; i<m_numBinsDoC; i++)
  {
    for (int j=0; j<m_numMolecPairs; j++)
    {
      if ( m_countIonsDoC[j] > 0)
      {
        m_DoCHist[i][3*j+1] /= m_DoCHist[i][3*j];
        m_DoCHist[i][3*j+2] /= m_DoCHist[i][3*j];
      }
      m_DoCHist[i][3*j] /= (numFrames * m_countIonsDoC[j] *  m_binSizeDoC);
    }
  }
  for (int i=0; i<m_numLayers; i++)
  {
    for (int j=0; j<m_numMolecPairs; j++)
    {
      m_coordNum[i][j] /= m_countCoordNum[i][j];
      for (int k=0; k<m_numBinsCoordNum; k++)
      {
        m_coordNumHist[i][k][j] /= m_countCoordNum[i][j];
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

/// Get address of first element in CoordNum
double* RDF::getCoordNumAddress(int i)
{
  return &(m_coordNum[i][0]);
}
/// Get address of first element in CoordNumHist
double* RDF::getCoordNumHistAddress(int i, int j)
{
  return &(m_coordNumHist[i][j][0]);
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
  const char * const headernames[] = { "z[A]",  "0",  "1",  "2",  "3",  "4",  "5",  "6",  "7",  "8",  "9",  "10",  "11",  "12",  "13", "14", "15", "16", "17", "18", "19", "20", "21", "22" };
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
  int varDim = 3*a_rdf->getNumMolecPairs();
  const char * const headernames[] = { "z[A]",  "DoC_0",  "cc_0", "cpc_0", "DoC_1",  "cc_1", "cpc_1",  "DoC_1",  "cc_2", "cpc_2",  "DoC_3",  "cc_3", "cpc_3",  "DoC_4",  "cc_4", "cpc_4",  "DoC_5",  "cc_5", "cpc_5",  "DoC_6",  "cc_6", "cpc_6",  "DoC_7",  "cc_7", "cpc_7",  "DoC_8",  "cc_8", "cpc_8",  "DoC_9",  "cc_9", "cpc_9",  "DoC_10",  "cc_10", "cpc_10",  "DoC_11",  "cc_11", "cpc_11",  "DoC_12",  "cc_12", "cpc_12",  "DoC_13",  "cc_13", "cpc_13" };
  double* data[500];
  for (int i=0; i<numBins; i++)
  {
    data[i] = a_rdf->getDoCHistAddress(i);
  }
  write_binned_data(a_filename, numBins, binSize, varDim, headernames, data);
  return a_filename;
}

const char* CoordNumHistWrite(RDF* a_rdf, const char* a_filename)
{
  double binSize = a_rdf->getBinSizeCoordNum();
  int numBins = a_rdf->getNumBinsCoordNum();
  int varDim = a_rdf->getNumMolecPairs();
  int numLayers = a_rdf->getNumLayers();
  const char * const headernames[] = { "nodeData" };
  double*** data = new double**[numLayers];
  for (int i=0; i<numLayers; i++)
  {
    data[i] = new double*[numBins];
    for (int j=0; j<numBins; j++)
    {
      data[i][j] =a_rdf->getCoordNumHistAddress(i,j);
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

const char* CoordNumWrite(RDF* a_rdf, const char* a_filename)
{
  int numLayers = a_rdf->getNumLayers();
  int varDim = a_rdf->getNumMolecPairs();
  double* layers;
  layers = new double[numLayers];
  a_rdf->getSystem().getLayerUpperBounds(numLayers,layers);
  const char * const headernames[] = { "z[A]", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10" };
  double** data;
  data = new double* [numLayers];
  for (int i=0; i<numLayers; i++)
  {
    data[i] = a_rdf->getCoordNumAddress(i);
  }
  write_layered_data(a_filename, numLayers, layers, varDim, headernames, data);
  delete data, layers;
  return a_filename;
}
