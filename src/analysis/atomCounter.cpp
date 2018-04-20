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
  m_binSize = 0.05;
  m_numBins = ceil(m_system.getBoxDim(2) / m_binSize );
  m_numLayers = m_system.getNumLayers();
  // resize vectors
  m_COMs.resize(m_system.getNumAtoms());
  m_numAtomsProfile.resize(m_numBins);
  m_numIonsProfile.resize(m_numBins);
  m_avgIonsInLayer.resize( m_numLayers );
  m_numAtomTypes=m_system.getNumAtomTypes();
  m_chargingParam.resize( m_system.getNumFrames() );
};

void AtomCounter::sample(const Frame& a_frame)
{
  int molecIndex = 0;
  int atomIndex = 0;
  int stepNum = a_frame.getStepNum();
  vector<array<int, NUM_ION_TYPES> > currentIonsInLayer;
  currentIonsInLayer.resize(m_numLayers);
  for (int i=0; i<m_system.getNumMolecTypes(); i++)
    {
      array<double , MAX_MEMBERS_PER_MOLEC > masses = m_system.getMassesOfType(i);
      for (int j=0; j < m_system.getNumMolecsOfType(i); j++)
	{
	  array<double, DIM > com;
	  int numMembers = m_system.getNumMembersMolec(i);
	  double totalMass = 0;
	  for (int k=0; k < numMembers; k++)
	    {
	      totalMass += masses[k];
	      // Get position of atom
	      array<double, DIM> position = a_frame.getAtom(atomIndex).getPosition();
	      // Bin atom by type
	      binAtom(position, i, k);
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
	  int* electrolyteID = new int;
	  if (m_system.isElectrolyte(i, electrolyteID))
	    {
	      binElectrolyteCOM(com, *electrolyteID);
	      countElectrolyteInLayer(com, currentIonsInLayer, *electrolyteID);
	    }
	  delete electrolyteID;
	  molecIndex++;
	}
    }
  for (int i=0; i<m_numLayers; i++)
    {
      for (int j=0; j<m_system.getNumElectrolyteSpecies(); j++)
	{
	  m_avgIonsInLayer[i][j] += currentIonsInLayer[i][j];
	}
    }
  //computeChargingParam(currentIonsInLayer);
}

void AtomCounter::binAtom(array<double, DIM>& a_position, int& a_molecType, int& a_molecMember)
{
  double pos_z = a_position[2];
  int bin = floor(pos_z / m_binSize);
  int atomType = m_system.getAtomType(a_molecType, a_molecMember);
#ifdef DEBUG
  if ( ! (atomType < m_system.getNumAtomTypes()) )
    {
      cout << a_molecType << " " << a_molecMember << " " << atomType << endl;
      cout << m_system.getNumAtomTypes() << endl;
    }
  assert ( atomType < m_system.getNumAtomTypes() ) ;
#endif
  m_numAtomsProfile[bin][atomType]++;
}

void AtomCounter::binElectrolyteCOM(array<double, DIM>& a_position, int& a_electrolyteID)
{
#ifdef DEBUG
  assert(a_electrolyteID > -1);
#endif
  double pos_z = a_position[2];
  int bin = floor(pos_z / m_binSize);
  m_numIonsProfile[bin][a_electrolyteID]++;
}

void AtomCounter::countElectrolyteInLayer(array<double, DIM>& a_position,  vector<array<int, NUM_ION_TYPES> >& a_currentIonsInLayer, int& a_electrolyteID)
{
#ifdef DEBUG
  assert(a_electrolyteID > -1);
#endif
  // Figure out what layer electrolyte molecule is in
  unsigned int layer = m_system.getLayer(a_position);
  // Increment count of molecule inside layer
  a_currentIonsInLayer[layer][a_electrolyteID]++;
}

const int AtomCounter::getNumAtomTypes()
{
    return m_numAtomTypes;
}

const int AtomCounter::getNumIonTypes()
{
    return NUM_ION_TYPES;
}

const int AtomCounter::getNumLayers()
{
    return m_numLayers;
}

void AtomCounter::normalize()
{
  // Normalize average ions in layer
  for (int i=0; i<m_numLayers; i++)
    {
      for (int j=0; j<m_system.getNumElectrolyteSpecies(); j++)
	{
	  m_avgIonsInLayer[i][j] /= m_system.getNumFrames();
	}
      for (int j=0; j<m_system.getNumAtomTypes(); j++)
	{
	  m_numAtomsProfile[i][j] /= m_system.getNumFrames();
	}
      for (int j=0; j<m_system.getNumMolecTypes(); j++)
	{
	  m_numIonsProfile[i][j] /= m_system.getNumFrames();
	}
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
}

const int AtomCounter::getNumBins() const
{
  return m_numBins;
}

const int AtomCounter::getBinSize() const
{
  return m_binSize;
}

double* AtomCounter::getACAtomsAddress(int i)
{
  return &(m_numAtomsProfile[i][0]);
}

double* AtomCounter::getACIonsAddress(int i)
{
  return &(m_numIonsProfile[i][0]);
}

double* AtomCounter::getACIonsLayersAddress(int i)
{
  return &(m_avgIonsInLayer[i][0]);
}

double AtomCounter::computeChargingParam(vector<array<int, NUM_ION_TYPES> >& a_ionsInLayer)
{

    // FIXME HERE!! COMPUTE CHARGING PARAMETER FOR CATHODE ANDA NODE SEPARATELY
    // FIXME        MAKE ROBUST WAY OF ACCESSING INDICES FOR CAT/ANODE AND CAT/ANIONS
    // FIXME        MICHELLE

  // N(V) total number of in-pore ions at charging voltage V
  int NV = a_ionsInLayer[0][0] + a_ionsInLayer[0][1];
  // N(V0) total number of in-pore ions at initial voltage V0
  int NV0 = 0.0;
  // Ncounter(V) number of in-pore counter-ions
  int NcounterV = a_ionsInLayer[0][0];
  int NcoV = a_ionsInLayer[0][1];
  // Nco(V) number of in-pore co-ions
  int NcounterV0 = 0.0;
  int NcoV0 = 0.0;
  double retVal = NcounterV-NcoV - NcounterV0+NcoV0;
  retVal = 1./retVal;
  retVal = retVal * (NV - NV0);
  return retVal;
}

const char* ACWriteDensity(AtomCounter* a_ac, const char* a_filename)
{
  double binSize = a_ac->getBinSize();
  int numBins = a_ac->getNumBins();
  int varDim = a_ac->getNumAtomTypes();
  const char * const headernames[] = { "nodeData" };
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

const char* ACWriteIons(AtomCounter* a_ac, const char* a_filename)
{
  double binSize = a_ac->getBinSize();
  int numBins = a_ac->getNumBins();
  int varDim = a_ac->getNumIonTypes();
  const char * const headernames[] = { "nodeData" };
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
  double binSize = a_ac->getBinSize();
  int numBins = a_ac->getNumLayers();
  int varDim = a_ac->getNumIonTypes();
  const char * const headernames[] = { "nodeData" };
  double** data;
  data = new double* [numBins];
  for (int i=0; i<numBins; i++)
    {
      data[i] = a_ac->getACIonsLayersAddress(i);
    }
  write_binned_data(a_filename, numBins, binSize, varDim, headernames, data);
  delete data;
  return a_filename;
}

const char* ACWriteCollectiveVars(AtomCounter* a_ac, const char* a_filename)
{
}
