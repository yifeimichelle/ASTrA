#include "frame.h"
#include "system.h"
#include "atom.h"
#include <vector>
#include "atomCounter.h"

AtomCounter::AtomCounter(System& a_system)
{
  m_system = a_system;
  m_numBins = 1000;
  // resize vectors
  m_COMs.resize(m_system.getNumAtoms());
  m_numAtomsProfile.resize(m_numBins);
  m_numIonsProfile.resize(m_numBins);
};

void AtomCounter::sample(Frame& a_frame)
{
  computeCOMs(a_frame);

  // ADD MORE STUFF HERE
  // FIXME MICHELLE
}

void AtomCounter::normalize()
{
}

void AtomCounter::print()
{
}

void AtomCounter::computeCOMs(Frame& a_frame)
{
  int indexCounter = 0;
  for (int i=0; i<m_system.getNumMolecTypes(); i++)
    {
      array<double , MAX_MEMBERS_PER_MOLEC > masses = m_system.getMassesOfType(i);
      for (int j=0; j<m_system.getNumMolecsOfType(i); j++)
	{
	  array<double, DIM > com;
	  for (int k=0; k<m_system.getNumMembersMolec(i); k++)
	    {
	      for (int l=0; l<DIM; l++)
		{
		  com[l] += masses[k];
		}
	    }
	  for (int l=0; l<DIM; l++)
	    {
	      com[l] /= m_system.getNumMembersMolec(i);
	    }
	  m_COMs[indexCounter]=com;
	  indexCounter++;
	}
    }
}

