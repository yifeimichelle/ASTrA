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
  m_binSize = 1;
  m_numBins = ceil(m_system.getBoxDim(2) / m_binSize );
  // resize vectors
  m_COMs.resize(m_system.getNumAtoms());
  m_numAtomsProfile.resize(m_numBins);
  m_numIonsProfile.resize(m_numBins);
  m_avgIonsInLayer.resize(m_system.getNumLayers() );
};

void AtomCounter::sample(const Frame& a_frame)
{
  int molecIndex = 0;
  int atomIndex = 0;
  vector<array<int, 3> > currentIonsInLayer;
  currentIonsInLayer.resize(m_system.getNumLayers());
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
  for (int i=0; i<m_system.getNumLayers(); i++)
    {
      for (int j=0; j<m_system.getNumElectrolyteSpecies(); j++)
	{
	  m_avgIonsInLayer[i][j] += currentIonsInLayer[i][j];
	}
    }
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

void AtomCounter::countElectrolyteInLayer(array<double, DIM>& a_position,  vector<array<int, 3> >& a_currentIonsInLayer, int& a_electrolyteID)
{
#ifdef DEBUG
  assert(a_electrolyteID > -1);
#endif
  // Figure out what layer electrolyte molecule is in
  unsigned int layer = m_system.getLayer(a_position);
  // Increment count of molecule inside layer
  a_currentIonsInLayer[layer][a_electrolyteID]++;
}

void AtomCounter::normalize()
{
  // Normalize average ions in layer
  for (int i=0; i<m_system.getNumLayers(); i++)
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
      for (int j=0; j<3; j++)
	{
	  cout << " " << m_numIonsProfile[i][j];
	}
      cout << endl;
    }
  cout << endl;
  for (int i=0; i<m_system.getNumLayers(); i++)
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
