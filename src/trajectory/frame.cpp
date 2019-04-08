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
  m_timestep = -1;
  m_numAtoms = m_system.getNumAtoms();
  m_numMolecules = m_system.getNumMolecules();
  m_every = m_system.getReadFrameEvery();

  // Initialize reading of trajectory file
  m_traj.open(m_system.getTrajFile().c_str());
  unsigned int trajNumAtoms;
  m_traj >> trajNumAtoms;
  assert(m_numAtoms==trajNumAtoms);
  m_atoms.resize(m_numAtoms);
  m_COMs.resize(m_numMolecules);
  m_traj.seekg(0, m_traj.beg);

  // Initialize frame charges from input charges
  setCharges(m_system);

  // Initialize reading of charges file
  if (m_system.hasChargeFile())
  {
    m_chg.open(m_system.getChargesFile().c_str());
    unsigned int chgNumAtoms;
    char tmp[256];
    m_chg >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp >> tmp;
    m_chg >> chgNumAtoms;
    m_numFluctuatingCharges = chgNumAtoms;
    m_chg.seekg(0, m_chg.beg);
  }

};
void Frame::skipStep()
{
  skipStep(m_every);
}
void Frame::skipStep(int a_every)
{
  // discard step from traj file
  for (int istep=0; istep < a_every-1; istep++)
  {
    char tmp[256];
    m_traj >> tmp >> tmp >> tmp >> tmp;
    for (unsigned int i=0; i<m_numAtoms; i++)
    {
      m_traj >> tmp >> tmp >> tmp >> tmp;
    }
  }
  m_totalStepNum++;
  char tmp[256];
  m_traj >> tmp >> tmp >> tmp >> tmp;
  string atomName;
  double x, y, z;
  for (unsigned int i=0; i<m_numAtoms; i++)
  {
    m_traj >> atomName >> x >> y >> z;
    m_atoms[i].setPosition(x,y,z);
    m_atoms[i].setName(atomName.c_str());
  }
}

void Frame::readZPStep()
{
  readZPStep(m_every);
}

void Frame::readZPStep(int a_every)
{
  int currentTimestep;
  // discard step from traj file
  for (int istep=0; istep < a_every-1; istep++)
  {
    char tmp[256];
    m_traj >> tmp >> tmp >> tmp >> tmp;
    for (unsigned int i=0; i<m_numAtoms; i++)
    {
      m_traj >> tmp >> tmp >> tmp >> tmp;
    }
  }
  // read step from traj file
  m_zpStepNum++;
  m_totalStepNum++;
  char tmp[256];
  m_traj >> tmp >> tmp >> tmp >> currentTimestep;
  m_timestep = currentTimestep;
  string atomName;
  double x, y, z;
  for (unsigned int i=0; i<m_numAtoms; i++)
  {
    m_traj >> atomName >> x >> y >> z;
    m_atoms[i].setPosition(x,y,z);
    m_atoms[i].setName(atomName.c_str());
  }
}

void Frame::readStep()
{
  readStep(m_every);
}

void Frame::readStep(int a_every)
{
  // discard step from traj file
  for (int istep=0; istep < a_every-1; istep++)
  {
    char tmp[256];
    m_traj >> tmp >> tmp >> tmp >> tmp;
    for (unsigned int i=0; i<m_numAtoms; i++)
    {
      m_traj >> tmp >> tmp >> tmp >> tmp;
    }
  }

  // read actual step
  m_stepNum++;
  m_totalStepNum++;
  char tmp[256];
  int currentTimestep;
  string atomName;
  double x, y, z;
  // read timestep
  m_traj >> tmp >> tmp >> tmp >> currentTimestep;
  // if current timestep is a duplicate, skip it
  // (ONLY IF READING EVERY STEP, otherwise throw error)
  //cout << currentTimestep << endl;
  while (currentTimestep <= m_timestep)
  {
    //cout << "current timestep is smaller than previous timestep" << endl;
    if (a_every > 1)
    {
      cerr << "ERROR: there is a duplicate or discontinuous timestep," << endl\
        << "and every " << a_every << " steps are being read." << endl \
        << "This is incompatible. m_timestep = " << m_timestep << ", " << endl\
        << "currentTimestep = " << currentTimestep << endl;
      exit(1);
      // or figure out how to handle this and break;
    }
    else
    {
      for (unsigned int i=0; i<m_numAtoms; i++)
      {
        m_traj >> atomName >> x >> y >> z;
        //m_atoms[i].setPosition(x,y,z); // don't need to read position
      }
      m_traj >> tmp >> tmp >> tmp >> currentTimestep;
    }
  }
  //cout << "reading next step" << endl;
  m_timestep = currentTimestep;
  for (unsigned int i=0; i<m_numAtoms; i++)
  {
    m_traj >> atomName >> x >> y >> z;
    m_atoms[i].setPosition(x,y,z);
    m_atoms[i].setName(atomName.c_str());
  }

}

void Frame::readCharges()
{
  char tmp[256];
  // skip lines
  for (unsigned int i=0; i<24; i++)
  {
    m_chg >> tmp;
  }
  int id;
  double q;
  for (unsigned int i=0; i<m_numFluctuatingCharges; i++)
  {
    m_chg >> id >> q;
    m_atoms[id-1].setCharge(q);
  }
}

void Frame::skipCharges()
{
  char tmp[256];
  // skip lines
  for (unsigned int i=0; i<24; i++)
  {
    m_chg >> tmp;
  }
  int id;
  double q;
  for (unsigned int i=0; i<m_numFluctuatingCharges; i++)
  {
    m_chg >> id >> q;
    m_atoms[id-1].setCharge(q);
  }
}

const int Frame::getTimestep() const
{
  return m_timestep;
}

const int Frame::getStepNum() const
{
  return m_stepNum;
}

const int Frame::getTotalStepNum() const
{
  return m_totalStepNum;
}

const int Frame::getZPStepNum() const
{
  return m_zpStepNum;
}

void Frame::clearFrame()
{
  for (unsigned int i=0; i<m_numAtoms; i++)
  {
    m_atoms[i].setPosition(-1, -1, -1);
    m_atoms[i].setName("");
  }
  for (unsigned int i=0; i<m_numMolecules; i++)
  {
    m_COMs[i].setPosition(-1, -1, -1);
    m_atoms[i].setName("");
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

const Atom& Frame::getAtomOfMolec(int a_molecIndex) const
{
  int atomIndex = m_system.getFirstAtomOfMolec(a_molecIndex);
  return m_atoms[atomIndex];
}


const double Frame::computeDistance(int a_i, int a_j) const
{
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

    retVal += dist*dist;
  }
  retVal = sqrt(retVal);
  return retVal;
}

const double Frame::computeMolecDistance(int a_i, int a_j) const
{
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
    retVal += dist*dist;
  }
  retVal = sqrt(retVal);
  // if (retVal == 0) //!!
  //   {
  //     cout << a_i << ", " << a_j << ": ";
  //     cout << posI[0] << " " ;
  //     cout << posI[1] << " " ;
  //     cout << posI[2] << " " << endl;
  //     cout << posJ[0] << " " ;
  //     cout << posJ[1] << " " ;
  //     cout << posJ[2] << " " << endl;
  //   }
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
    char str[100];
    sprintf(str, "%d", i);
    m_COMs[i].setName(str);
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
const unsigned int Frame::getCurrentNumAtomsInLayer(int a_layerIdx, int a_molID) const
{
  return m_atomLayers[a_layerIdx][a_molID].size();
}

const unsigned int Frame::getCurrentNumMolecsInLayer(int a_layerIdx, int a_molID) const
{
  return m_COMLayers[a_layerIdx][a_molID].size();
}

void Frame::setCharges(const System& a_system)
{
  unsigned int atomIdx = 0;
  double q;
  for (unsigned int i=0; i<a_system.getNumMolecTypes(); i++)
  {
    array<double , MAX_MEMBERS_PER_MOLEC > charges = a_system.getChargesOfType(i);
    unsigned int numMolecs=a_system.getNumMolecsOfType(i);
    for (unsigned int j=0; j<numMolecs; j++)
    {
      unsigned int numMembers=a_system.getNumMembersMolec(i);
      for (unsigned int k=0; k<numMembers; k++)
      {
        m_atoms[atomIdx].setCharge(charges[k]);
        atomIdx++;
      }
    }
  }
}

const double Frame::sumCharges(const System& a_system, int a_molID) const
{
  double sum = 0;
  unsigned int atomIdx = 0;
  for (unsigned int i=0; i<a_system.getNumMolecTypes(); i++)
  {
    unsigned int numMolecs=a_system.getNumMolecsOfType(i);
    for (unsigned int j=0; j<numMolecs; j++)
    {
      unsigned int numMembers=a_system.getNumMembersMolec(i);
      for (unsigned int k=0; k<numMembers; k++)
      {
        if ( i == a_molID-1 )
        {
          sum += m_atoms[atomIdx].getCharge();
        }
        atomIdx++;
      }
    }
  }
  return sum;
}
