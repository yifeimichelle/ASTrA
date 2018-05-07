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
  m_totalStepNum = -1;
  m_zpStepNum = -1;
  m_numAtoms = a_system.getNumAtoms();
  m_numMolecules = a_system.getNumMolecules();
  m_traj.open(a_system.getTrajFile());
  unsigned int trajNumAtoms;
  m_traj >> trajNumAtoms;
  assert(m_numAtoms==trajNumAtoms);
  m_atoms.resize(m_numAtoms);
  m_COMs.resize(m_numMolecules);
  m_traj.seekg(0, m_traj.beg);
};

void Frame::skipStep()
{
  m_totalStepNum++;
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

void Frame::readZPStep()
{
  // read step from traj file
  m_zpStepNum++;
  m_totalStepNum++;
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

void Frame::readStep()
{
  // read step from traj file
  m_stepNum++;
  m_totalStepNum++;
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

const unsigned int Frame::getTotalStepNum() const
{
  return m_totalStepNum;
}

const unsigned int Frame::getZPStepNum() const
{
  return m_zpStepNum;
}

void Frame::clearFrame()
{
  for (unsigned int i=0; i<m_numAtoms; i++)
    {
      m_atoms[i].setPosition(-1, -1, -1);
    }
    for (unsigned int i=0; i<m_numMolecules; i++)
    {
      m_COMs[i].setPosition(-1, -1, -1);
    }
  for (unsigned int i=0; i<NUM_LAYERS; i++)
    {
      for (unsigned int j=0; j<MAX_NUM_TYPES; j++)
	{
	  m_atomLayers[i][j].clear();
	  m_COMLayers[i][j].clear();
	  m_ZPatomLayers[i][j].clear();
	  m_ZPCOMLayers[i][j].clear();
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

const Atom& Frame::getMolec(int a_molecIndex) const
{
  return m_COMs[a_molecIndex];
}

Atom& Frame::getMolec(int a_molecIndex)
{
  return m_COMs[a_molecIndex];
} 

const double Frame::computeDistance(int a_i, int a_j) const
{
  //array<double, DIM > boxDims = m_system.getBoxDims();
  array<double, DIM > posI = m_atoms[a_i].getPosition();
  array<double, DIM > posJ = m_atoms[a_j].getPosition();
  double retVal = 0.0;
  double dim;
  for (int i=0; i<DIM; i++)
    {
      dim = m_system.getBoxDim(i);
      double dist = posI[i] - posJ[i];
      if(m_system.isPeriodic(i))
	{
	  dist -= round(dist/dim) * dim;
	}
      // if(m_system.isPeriodic(i))
      // 	{
      // 	  dist -= boxDims[i]*floor(dist/(0.5*boxDims[i]));
      // 	}
      retVal += dist*dist;
    }
  retVal = sqrt(retVal);
  return retVal;
}

const double Frame::computeMolecDistance(int a_i, int a_j) const
{
  //array<double, DIM > boxDims = m_system.getBoxDims();
  array<double, DIM > posI = m_COMs[a_i].getPosition();
  array<double, DIM > posJ = m_COMs[a_j].getPosition();
  double retVal = 0.0;
  double dim;
  for (int i=0; i<DIM; i++)
    {
      dim = m_system.getBoxDim(i);
      double dist = posI[i] - posJ[i];
      if(m_system.isPeriodic(i))
	{
	  dist -= round(dist/dim) * dim;
	}
      // if(m_system.isPeriodic(i))
      // 	{
      // 	  dist -= boxDims[i]*floor(dist/(0.5*boxDims[i]));
      // 	}
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
  for (int i=0; i<m_numMolecules; i++)
    {
      m_COMs[i].setPosition(a_COMs[i][0],a_COMs[i][1],a_COMs[i][2]);
    }
}


void Frame::assignZPAtomToLayer(unsigned int a_index, unsigned int a_type, unsigned int a_layer)
{
  m_ZPatomLayers[a_layer][a_type].push_back(a_index);
}

void Frame::assignZPIonToLayer(unsigned int a_index, unsigned int a_type, unsigned int a_layer)
{
  m_ZPCOMLayers[a_layer][a_type].push_back(a_index);
}

void Frame::assignAtomToLayer(unsigned int a_index, unsigned int a_type, unsigned int a_layer)
{
  m_atomLayers[a_layer][a_type].push_back(a_index);
}

void Frame::assignIonToLayer(unsigned int a_index, unsigned int a_type, unsigned int a_layer)
{
  m_COMLayers[a_layer][a_type].push_back(a_index);
}

vector<int>* Frame::getAtomsInLayer(int a_layerIdx) const
{
  return (vector<int>* )&m_atomLayers[a_layerIdx][0];
}

vector<int>* Frame::getMoleculesInLayer(int a_layerIdx) const
{
  return (vector<int>* )&m_COMLayers[a_layerIdx][0];
}


void Frame::printAtomsInLayer(unsigned int a_layer)
{
  cout << "Atom indices in layer " << a_layer <<":" << endl;
  for (int type = 0; type < m_system.getNumAtomTypes(); type++)
    {
      cout << "  Type " << type << ":";
      for (vector<int>::iterator it = m_atomLayers[a_layer][type].begin(); it != m_atomLayers[a_layer][type].end(); it++)
	{
	  cout << " " << *it;
	}
      cout << endl;
    }
}

void Frame::printAtomsInLayerCheck(unsigned int a_layer)
{
  array<vector<int >, MAX_NUM_TYPES > retVals;
  for (int type=0; type < m_system.getNumAtomTypes(); type++)
    {
    for (int atom=0; atom<m_system.getNumOfType(type); atom++)
      {
	int index = m_system.getIndexOfType(type, atom);
	int layer = getLayerOf(index);
	if (layer == a_layer)
	  {
	    retVals[type].push_back(index);
	  }
      }
    }
  cout << "Atom indices in layer " << a_layer << " (debug purposes):" << endl;
  for (int type = 0; type < m_system.getNumAtomTypes(); type++)
    {
      cout << "  Type " << type << ":";
      for (vector<int>::iterator it = retVals[type].begin(); it != retVals[type].end(); it++)
	{
	  cout << " " << *it;
	}
      cout << endl;
    }

}

void Frame::printMolecsInLayer(unsigned int a_layer)
{
  cout << "Molecule indices in layer " << a_layer <<":" << endl;
  for (int type = 0; type < m_system.getNumMolecTypes(); type++)
    {
      cout << "  Type " << type << ":";
      for (vector<int>::iterator it = m_COMLayers[a_layer][type].begin(); it != m_COMLayers[a_layer][type].end(); it++)
	{
	  cout << " " << *it;
	}
      cout << endl;
    }
}



