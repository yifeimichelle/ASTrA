#include "frame.h"
#include "system.h"
#include "atom.h"
#include <vector>
#include "atomCounter.h"
#include <assert.h>
#include <iostream>
#include "writer.h"
 
using namespace std;

AtomCounter::AtomCounter(System& a_system)
{
  m_system = a_system;
  m_saveFrameEvery = 1;
  m_binSize = 0.05;
  m_numBins = ceil(m_system.getBoxDim(2) / m_binSize );
  m_numLayers = m_system.getNumLayers();
  // resize vectors
  m_COMs.resize(m_system.getNumMolecules());
  m_numAtomsProfile.resize(m_numBins);
  m_densityProfile.resize(m_numBins);
  m_numIonsProfile.resize(m_numBins);
  m_avgIonsInLayer.resize( m_numLayers );
  m_numSavedFrames = ceil(m_system.getNumFrames() / m_saveFrameEvery) + 1;
  m_numIonsInLayerTime.resize(m_numSavedFrames);
  for (int i=0; i<m_numSavedFrames; i++)
    {
      m_numIonsInLayerTime[i].resize( m_numLayers );
    }
  m_numAtomTypes=m_system.getNumAtomTypes();
  m_chargingParam.resize( m_system.getNumFrames() );
};

void AtomCounter::sample(Frame& a_frame)
{
  int molecIndex = 0;
  int atomIndex = 0;
  int isElectrolyte = 0;
  array<double, DIM > com;
  int stepNum = a_frame.getStepNum();
  vector<array<int, NUM_ION_TYPES> > currentIonsInLayer;
  currentIonsInLayer.resize(m_numLayers);
  for (int i=0; i<m_system.getNumMolecTypes(); i++)
    {
      int numMembers = m_system.getNumMembersMolec(i);
      int* electrolyteID = new int;
      if (m_system.isElectrolyte(i, electrolyteID))
	{
	  isElectrolyte = 1;
	}
      else
	{
	  isElectrolyte = 0;
	}

      array<double , MAX_MEMBERS_PER_MOLEC > masses = m_system.getMassesOfType(i);
      double totalMass = 0;
      for (int k=0; k < numMembers; k++)
	{
	  totalMass += masses[k];
	}
      for (int j=0; j < m_system.getNumMolecsOfType(i); j++)
	{
	  com.fill(0); // fill with zeros, otherwise will keep same data as before
	  for (int k=0; k < numMembers; k++)
	    {
	      // Get position of atom
	      array<double, DIM> position = a_frame.getAtom(atomIndex).getPosition();
	      // Bin atom by type and add to density
	      binAtom(a_frame, atomIndex, position, i, k, masses[k], isElectrolyte);
	      // Compute center of mass
	      for (int l=0; l < DIM; l++)
		{
		  com[l] += position[l]*masses[k];
		}
	      atomIndex++;
	    }
	  for (int l=0; l<DIM; l++)
	    {
	      com[l] /= totalMass;
	    }
	  m_COMs[molecIndex]=com;
	  binElectrolyteCOM(a_frame, molecIndex, com, i, currentIonsInLayer, *electrolyteID, isElectrolyte);
	  molecIndex++;
	}
      delete electrolyteID;
    }
  for (int i=0; i<m_numLayers; i++)
    {
      for (int j=0; j<m_system.getNumElectrolyteSpecies(); j++)
	{
	  m_avgIonsInLayer[i][j] += currentIonsInLayer[i][j];
	}
    }
  if ( (a_frame.getStepNum() % m_saveFrameEvery == 0) || ( a_frame.getStepNum() == m_system.getNumFrames() ) )
    {
      for (int i=0; i<m_numLayers; i++)
	{
	  for (int j=0; j<m_system.getNumElectrolyteSpecies(); j++)
	    {
	      int timeIndex = (ceil)( ((double)a_frame.getStepNum()) / ((double)m_saveFrameEvery) );
	      m_numIonsInLayerTime[timeIndex][i][j] = currentIonsInLayer[i][j];
	    }
	}
    }
  a_frame.setCOMs(m_COMs);
  // compute charging parameter
}

void AtomCounter::binAtom(Frame& a_frame,  int& a_atomIndex, array<double, DIM>& a_position, int& a_molecType, int& a_molecMember, double& a_mass, int& a_isElectrolyte)
{
  double pos_z = a_position[2];
  int bin = floor(pos_z / m_binSize);
  int atomType = m_system.getAtomType(a_molecType, a_molecMember);
#ifdef DEBUG
  assert ( atomType < m_system.getNumAtomTypes() ) ;
#endif
  // Increment number of atoms in bin
  m_numAtomsProfile[bin][atomType]++;
  if (a_isElectrolyte)
    {
      // Increment density of electrolyte in bin
      m_densityProfile[bin] += a_mass;
    }
  unsigned int layer = m_system.getLayer(a_position);
  a_frame.assignAtomToLayer(a_atomIndex, atomType, layer);
}

void AtomCounter::binElectrolyteCOM(Frame& a_frame, int& a_molecIndex, array<double, DIM>& a_position, int& a_molecType, vector<array<int, NUM_ION_TYPES> >& a_currentIonsInLayer, int& a_electrolyteID, int& a_isElectrolyte)
{
  unsigned int layer = m_system.getLayer(a_position);
  a_frame.assignAtomToLayer(a_molecIndex, a_molecType, layer);
  if(a_isElectrolyte)
    {
#ifdef DEBUG
      assert(a_electrolyteID > -1);
      assert(a_electrolyteID < 3);
#endif
      double pos_z = a_position[2];
      int bin = floor(pos_z / m_binSize);
#ifdef DEBUG
      cout << "bin " << bin << " of " << m_numBins << " bins" << endl;
      double x = pos_z / m_system.getBoxDim(2);
      cout << "z position " << pos_z << " of electrolyte " << a_electrolyteID << endl;
      assert( bin < m_numBins );
#endif
      m_numIonsProfile[bin][a_electrolyteID]++;
      // Figure out what layer electrolyte molecule is in
      // Increment count of molecule inside layer
      a_currentIonsInLayer[layer][a_electrolyteID]++;
    }
}

const int AtomCounter::getNumAtomTypes()
{
    return m_numAtomTypes;
}

const int AtomCounter::getNumIonTypes()
{
    return NUM_ION_TYPES;
}

const int AtomCounter::getNumLayers() const
{
    return m_numLayers;
}

void AtomCounter::normalize()
{
  // Normalize average ions in layer:
  double avogadro = 6.022140857;
  // Divide by number of frames read
  double numFrames = m_system.getNumFrames();
  // Convert density to g/ml
  double normDensity = numFrames * m_binSize * m_system.getBoxDim(0) * m_system.getBoxDim(1) * 10. / avogadro;
  for (int i=0; i<m_numLayers; i++)
    {
      for (int j=0; j<m_system.getNumElectrolyteSpecies(); j++)
	{
	  m_avgIonsInLayer[i][j] /= numFrames;
	}
    }
    for (int i=0; i<m_numBins; i++)
      {
	for (int j=0; j<m_system.getNumAtomTypes(); j++)
	  {
	    m_numAtomsProfile[i][j] /= numFrames;
	  }
	for (int j=0; j<getNumIonTypes(); j++)
	  {
	    m_numIonsProfile[i][j] /= numFrames;
	  }
	m_densityProfile[i] /= normDensity ;
      }
}


void AtomCounter::print()
{
  for (int i=0; i<m_numBins; i++)
    {
      cout << i*getBinSize();
      for (int j=0; j<m_system.getNumAtomTypes(); j++)
	{
	  cout << " " << m_numAtomsProfile[i][j];
	}
      for (int j=0; j<NUM_ION_TYPES; j++)
	{
	  cout << " " << m_numIonsProfile[i][j];
	}
      cout << endl;
    }
  cout << endl;
  for (int i=0; i<m_numLayers; i++)
    {
      for (int j=0; j<m_system.getNumElectrolyteSpecies(); j++)
	{
	  cout << m_avgIonsInLayer[i][j] << " ";
	}
      cout << endl;
    }
  cout << endl;
  for (int i=0; i<m_numSavedFrames; i++)
    {
      for (int j=0; j<m_numLayers; j++)
  	{
  	  for (int k=0; k<m_system.getNumElectrolyteSpecies(); k++)
  	    {
  	      cout << m_numIonsInLayerTime[i][j][k] << " ";
  	    }
  	  cout << endl;
  	}
    }
}

void AtomCounter::printDensity()
{
  for (int i=0; i<m_numBins; i++)
    {
      cout << i*getBinSize() << " " << m_densityProfile[i] << endl;
    }
  cout << endl;
}

const int AtomCounter::getNumBins() const
{
  return m_numBins;
}

const double AtomCounter::getBinSize() const
{
  return m_binSize;
}

double* AtomCounter::getACAtomsAddress(int i)
{
  return &(m_numAtomsProfile[i][0]);
}

double* AtomCounter::getACDensityAddress(int i)
{
  return &(m_densityProfile[i]);
}

double* AtomCounter::getACIonsAddress(int i)
{
  return &(m_numIonsProfile[i][0]);
}

double* AtomCounter::getACIonsLayersAddress(int i)
{
  return &(m_avgIonsInLayer[i][0]);
}

double* AtomCounter::getACIonsLayersTimeAddress(int i, int j)
{
  return &(m_numIonsInLayerTime[i][j][0]);
}


double AtomCounter::computeChargingParam(vector<array<int, NUM_ION_TYPES> >& a_ionsInLayer)
{

    // FIXME HERE!! COMPUTE CHARGING PARAMETER FOR CATHODE AND ANODE SEPARATELY
    // FIXME        MAKE ROBUST WAY OF ACCESSING INDICES FOR CAT/ANODE AND CAT/ANIONS
    // FIXME        MICHELLE

  // N(V) total number of in-pore ions at charging voltage V
  int NV = a_ionsInLayer[0][0] + a_ionsInLayer[0][1];
  // N(V0) total number of in-pore ions at initial voltage V0
  int NV0 = 0.0;
  // Nco,Ncounter(V) number of in-pore co- and counter-ions at charging voltage V
  int NcounterV = a_ionsInLayer[0][0];
  int NcoV = a_ionsInLayer[0][1];
  // Nco,Ncounter(V0) number of in-pore co- and counter-ions at initial voltage V0
  int NcounterV0 = 0.0;
  int NcoV0 = 0.0;
  double retVal = NcounterV-NcoV - NcounterV0+NcoV0;
  retVal = 1./retVal;
  retVal = retVal * (NV - NV0);
  return retVal;
}

const System& AtomCounter::getSystem() const
{
  return m_system;
}

const int AtomCounter::getSaveFrameInterval() const
{
  return m_saveFrameEvery;
}

const int AtomCounter::getNumSavedFrames() const
{
  return m_numSavedFrames;
}

// const vector<array<double, DIM >  > AtomCounter::getCOMs() const
// {
//   return m_COMs;
// }


const char* ACWriteAtomCounts(AtomCounter* a_ac, const char* a_filename)
{
  double binSize = a_ac->getBinSize();
  int numBins = a_ac->getNumBins();
  int varDim = a_ac->getNumAtomTypes();
  const char * const headernames[] = { "z[A]",  "0",  "1",  "2",  "3",  "4",  "5",  "6",  "7",  "8",  "9",  "10",  "11",  "12",  "13", "14" };
  double** data;
  data = new double* [numBins];
  for (int i=0; i<numBins; i++)
    {
      data[i] = a_ac->getACAtomsAddress(i);
    }
  write_binned_data(a_filename, numBins, binSize, varDim, headernames, data);
  delete data;
  return a_filename;
}

const char* ACWriteDensity(AtomCounter* a_ac, const char* a_filename)
{
  double binSize = a_ac->getBinSize();
  int numBins = a_ac->getNumBins();
  int varDim = 1;
  const char * const headernames[] = { "z[A]", "dens[g/ml]" };
  double** data;
  data = new double* [numBins];
  for (int i=0; i<numBins; i++)
    {
      data[i] = a_ac->getACDensityAddress(i);
    }
  write_binned_data(a_filename, numBins, binSize, varDim, headernames, data);
  delete data;
  return a_filename;
}

const char* ACWriteIons(AtomCounter* a_ac, const char* a_filename)
{
  double binSize = a_ac->getBinSize();
  int numBins = a_ac->getNumBins();
  int varDim = a_ac->getNumIonTypes();
  const char * const headernames[] = { "z[A]", "cation", "anion", "solvent" };
  double** data;
  data = new double* [numBins];
  for (int i=0; i<numBins; i++)
    {
      data[i] = a_ac->getACIonsAddress(i);
    }
  write_binned_data(a_filename, numBins, binSize, varDim, headernames, data);
  delete data;
  return a_filename;
}

const char* ACWriteIonsInLayers(AtomCounter* a_ac, const char* a_filename)
{
  int numLayers = a_ac->getNumLayers();
  double* layers;
  layers = new double[numLayers];
  a_ac->getSystem().getLayerUpperBounds(numLayers,layers);
  int varDim = a_ac->getNumIonTypes();
  const char * const headernames[] = { "z[A]", "cation", "anion", "solvent" };
  double** data;
  data = new double* [numLayers];
  for (int i=0; i<numLayers; i++)
    {
      data[i] = a_ac->getACIonsLayersAddress(i);
    }
  write_layered_data(a_filename, numLayers, layers, varDim, headernames, data);
  delete data, layers;
  return a_filename;
}

const char* ACWriteIonsInLayersTime(AtomCounter* a_ac, const char* a_filename)
{
  int numFrames = a_ac->getNumSavedFrames();
  int saveFrameEvery = a_ac->getSaveFrameInterval();
  int numLayers = a_ac->getNumLayers();
  int varDim = a_ac->getNumIonTypes();
  const char * const headernames[] = { "t", "data" };
  double*** data;
  data = new double** [numFrames];
  for (int i=0; i<numFrames; i++)
    {
      // FIXME need better way to access pointers than doing this every time, it is very tedious
      data[i] = new double* [numLayers];
      for (int j=0; j<numLayers; j++)
	{
	  data[i][j] = a_ac->getACIonsLayersTimeAddress(i,j);
	}
    }
  write_layered_time_data_by_var(a_filename, numFrames, saveFrameEvery, numLayers, varDim, headernames, data);
  delete data;
  return a_filename;
}

const char* ACWriteCollectiveVars(AtomCounter* a_ac, const char* a_filename)
{
}
