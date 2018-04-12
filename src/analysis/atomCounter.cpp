#include "frame.h"
#include "system.h"
#include "atom.h"
#include <vector>
#include "atomCounter.h"
#include <assert.h>
#include <iostream>

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
};

void AtomCounter::sample(Frame& a_frame)
{
  computeCOMs(a_frame);
  computeDensity(a_frame);
  // ADD MORE STUFF HERE
  // FIXME MICHELLE
  int molecIndex = 0;
  int atomIndex = 0;
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
	      binAtom(position, i, j);
	      // Compute center of mass
	      for (int l=0; l < DIM; l++)
		{
		  com[l] += position[l]*masses[k];
		}
	      atomIndex++;
	    }
	  for (int l=0; l<DIM; l++)
	    {
	      com[l] /= (numMembers*totalMass);
	    }
	  m_COMs[molecIndex]=com;
	  int* electrolyteID = new int;
	  if (m_system.isElectrolyte(i, electrolyteID))
	    {
	      binElectrolyteCOM(com, *electrolyteID);
	    }
	  delete electrolyteID;
	  molecIndex++;
	}
    }
}

void AtomCounter::binAtom(array<double, DIM> a_position, int a_molecType, int a_molecMember)
{
  double pos_z = a_position[2];
  int bin = floor(pos_z / m_binSize);
  int atomType = m_system.getAtomType(a_molecType, a_molecMember);
  m_numAtomsProfile[bin][atomType]++;
}

void AtomCounter::binElectrolyteCOM(array<double, DIM> a_position, int a_electrolyteID)
{
  assert(a_electrolyteID > -1);
  double pos_z = a_position[2];
  int bin = floor(pos_z / m_binSize);
  m_numIonsProfile[bin][a_electrolyteID]++;
}

void AtomCounter::normalize()
{
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
}

const int AtomCounter::getNumBins() const
{
  return m_numBins;
}

const int AtomCounter::getBinSize() const
{
  return m_binSize;
}


void AtomCounter::computeCOMs(Frame& a_frame)
{
}

void AtomCounter::computeDensity(Frame& a_frame)
{ 
}
