#include <fstream>
#include <istream>
#include <iostream>
#include <cassert>
#include "math.h"
#include "atom.h"
#include "system.h"
#include "frame.h"

using namespace std;

Frame::Frame() {
};
Frame::Frame(System& a_system)
{
  m_system = a_system;
  m_stepNum = -1;
  m_numAtoms = a_system.getNumAtoms();
  m_numMolecs = a_system.getNumMolecules();
  m_traj.open(a_system.getTrajFile());
  unsigned int trajNumAtoms;
  m_traj >> trajNumAtoms;
  assert(m_numAtoms==trajNumAtoms);
  m_atoms.resize(m_numAtoms);
  m_COMs.resize(m_numMolecs);
  m_traj.seekg(0, m_traj.beg);
};
void Frame::readStep()
{
  // read step from traj file
  m_stepNum++;
  char tmp[256];
  m_traj >> tmp >> tmp >> tmp >> tmp;
  string atomName;
  double x, y, z;
  for (unsigned int i=0; i<m_numAtoms; i++)
    {
      m_traj >> atomName >> x >> y >> z;
      m_atoms[i].setPosition(x,y,z);
    }
}

const unsigned int Frame::getStepNum() const
{
  return m_stepNum;
}

void Frame::clearFrame()
{
  for (unsigned int i=0; i<m_numAtoms; i++)
    {
      m_atoms[i].setPosition(-1, -1, -1);
    }
    for (unsigned int i=0; i<m_numMolecs; i++)
    {
      m_COMs[i].setPosition(-1, -1, -1);
    }
  for (unsigned int i=0; i<NUM_LAYERS; i++)
    {
      for (unsigned int j=0; j<MAX_NUM_TYPES; j++)
	{
	  m_atomLayers[i][j].clear();
	  m_COMLayers[i][j].clear();
	}
    }
}

const Atom& Frame::getAtom(int a_atomIndex) const
{
  return m_atoms[a_atomIndex];
}

Atom& Frame::getAtom(int a_atomIndex)
{
  return m_atoms[a_atomIndex];
} 

const double Frame::computeDistance(int a_i, int a_j) const
{
  array<double, DIM > boxDims = m_system.getBoxDims();
  array<double, DIM > posI = m_atoms[a_i].getPosition();
  array<double, DIM > posJ = m_atoms[a_j].getPosition();
  double retVal = 0.0;
  for (int i=0; i<DIM; i++)
    {
      double dist = abs(posI[i] - posJ[i]);
      if(m_system.isPeriodic(i))
	{
	  dist -= boxDims[i]*floor(dist/(0.5*boxDims[i]));
	}
      retVal += dist*dist;
    }
  retVal = sqrt(retVal);
  return retVal;
}

const double Frame::computeMolecDistance(int a_i, int a_j) const
{
  array<double, DIM > boxDims = m_system.getBoxDims();
  array<double, DIM > posI = m_COMs[a_i].getPosition();
  array<double, DIM > posJ = m_COMs[a_j].getPosition();
  double retVal = 0.0;
  for (int i=0; i<DIM; i++)
    {
      double dist = abs(posI[i] - posJ[i]);
      if(m_system.isPeriodic(i))
	{
	  dist -= boxDims[i]*floor(dist/(0.5*boxDims[i]));
	}
      retVal += dist*dist;
    }
  retVal = sqrt(retVal);
  return retVal;
}

const unsigned int Frame::getLayerOf(unsigned int a_index) const
{
  array<double, DIM > position = m_atoms[a_index].getPosition();
  unsigned int retval = m_system.getLayer(position);
  return retval;
}

const unsigned int Frame::getLayerOfMolec(unsigned int a_index) const
{
  array<double, DIM > position = m_COMs[a_index].getPosition();
  unsigned int retval = m_system.getLayer(position);
  return retval;
}

void Frame::setCOMs(vector<array<double, DIM >  > a_COMs)
{
  for (int i=0; i<m_numMolecs; i++)
    {
      m_COMs[i].setPosition(a_COMs[i][0],a_COMs[i][1],a_COMs[i][2]);
    }
}


void Frame::assignAtomToLayer(unsigned int a_index, unsigned int a_type, unsigned int a_layer)
{
  m_atomLayers[a_layer][a_type].push_back(a_index);
}

void Frame::assignIonToLayer(unsigned int a_index, unsigned int a_type, unsigned int a_layer)
{
  m_COMLayers[a_layer][a_type].push_back(a_index);
}
