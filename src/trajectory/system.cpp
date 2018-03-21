#include <math.h>
#include <iostream>
#include <cstring>
#include <array>
#include <fstream>
#include <string>
#include "system.h"

using namespace std;

unsigned int System::getAtomType(unsigned int a_molecType, unsigned int a_molecAtom)
{
  unsigned int retVal = 0;
  for (int i=0; i<a_molecType; i++)
    {
      for (int j=0; j<m_numAtomsMolec[i]; j++)
	{
	  retVal++;
	}
    }
  for (int j=0; j<a_molecAtom; j++)
    {
      retVal++;
    }
  return retVal;
}

System::System()
{
};

System::System(const string& a_inputFile)
{
    ifstream system(a_inputFile.c_str());
    cout << "Reading input file ";
    cout << a_inputFile << endl;

    system >> m_trajFile;
    system >> m_numFrames;

    for (int i=0; i<DIM; i++) {
        system >> m_boxDims[i];
    }
    for (int i=0; i<DIM; i++) {
        system >> m_boxPeriodic[i];
    }
    cout << m_trajFile << endl;
    cout << "Box dims: ";
    for (int i=0; i<DIM; i++) {
        cout << m_boxDims[i] << " ";
    }
    cout << endl;
    system >> m_lowerElecTop >> m_upperElecBot;
    cout << "Electrode limits: " << m_lowerElecTop << " " << m_upperElecBot << endl;
    system >> m_numMolecTypes;
    m_numAtomTypes = 0;
    for (unsigned int i=0; i<m_numMolecTypes; i++) {
        system >> m_numMolecs[i] >> m_numAtomsMolec[i];
	m_numAtomTypes += m_numAtomsMolec[i];
    }
    system >> m_cationID >> m_anionID >> m_lowerElecID >> m_upperElecID ;
    system >> m_anodeIsLower ;
    if (m_anodeIsLower)
      {
	m_anodeID = m_lowerElecID;
	m_cathodeID = m_upperElecID;
      }
    else
      {
	m_anodeID = m_upperElecID;
	m_cathodeID = m_lowerElecID;
      }
    system >> m_solventID ;
    system >> m_capID ;
    if (m_solventID == 0) {
      m_boolWithSolvent = 0; }
    else {
      m_boolWithSolvent = 1; }

    if (m_capID == 0) {
      m_boolWithCap = 0; }
    else {
      m_boolWithCap = 0; }
    cout << "number of atom types: " << m_numAtomTypes << endl;
    system >> m_stepInterval >> m_stepTime;
    system >> m_numPairs;
    m_rdfPairs.resize(m_numPairs);
    cout << "Pair correlations:" << endl;
    for (unsigned int i=0; i<m_numPairs; i++) {
      unsigned int molecA, atomA, molecB, atomB;
      system >> molecA >> atomA >> molecB >> atomB;
      unsigned int atomTypeA = getAtomType(molecA-1, atomA-1);
      unsigned int atomTypeB = getAtomType(molecB-1, atomB-1);
      cout << i << " " << atomTypeA << " " << atomTypeB << endl;
      m_rdfPairs[i] = make_pair(atomTypeA, atomTypeB);
    }

    m_frameTime = m_stepInterval * m_stepTime;
    m_typeAtomIndices.resize(m_numAtomTypes);
    unsigned int counter=0;
    for (int i=0; i<m_numMolecTypes; i++)
      {
	for (int j=0; j<m_numMolecs[i]; j++)
	  {
	    for (int k=0; k<m_numAtomsMolec[i]; k++)
	      {
		m_typeAtomIndices[getAtomType(i,k)].push_back(counter);
		counter++;
	      }
	  }
      }
    
};

const int System::getNumAtoms() const
{
    int retVal = 0;
    for (unsigned int i=0; i<m_numMolecTypes; i++) {
        retVal += m_numMolecs[i]*m_numAtomsMolec[i];
    }
    return retVal;
}

const float System::getFrameTime() const
{
    return m_frameTime;
}

const string& System::getTrajFile() const
{
    return m_trajFile;
}

void System::printPairCorrelations() const
{
  for (int i=0; i<m_numPairs; i++)
    {
      cout << "pair " << i << " : ";
      cout << m_rdfPairs[i].first << " " << m_rdfPairs[i].second << endl;
      
    }
}

void System::printTypeAtomIndices() const
{
  for (int i=0; i<m_numAtomTypes; i++)
    {
      cout << "type " << i << " : " << endl;
      for (int j=0; j<m_typeAtomIndices[i].size() ; j++)
	{
	  cout << m_typeAtomIndices[i][j] << " ";
	}
      cout << endl;
    }
}

const unsigned int System::getNumOfType(unsigned int a_type) const
{
  return m_typeAtomIndices[a_type].size();
}

const unsigned int System::getIndexOfType(unsigned int a_type, unsigned int a_idx) const
{
  return m_typeAtomIndices[a_type][a_idx];
}

const unsigned int System::getNumPairs() const
{
  return m_numPairs;
}

const pair<unsigned int, unsigned int > System::getPairCorrelation(unsigned int a_pair) const
{
  return m_rdfPairs[a_pair];
}

const unsigned int System::getNumFrames() const
{
  return m_numFrames;
}

const array<double, DIM > System::getBoxDims() const
{
  return m_boxDims;
}
const unsigned int System::isPeriodic(int i) const
{
  return m_boxPeriodic[i];
}

const unsigned int System::getNumLayers() const
{
  return 3;
}

const unsigned int System::getLayer(array<double, DIM>& a_position) const
{
  unsigned int retVal;
  double z=a_position[DIM-1];
  if (z < m_lowerElecTop )
    {
      retVal = 0;
    }
  else if (z < m_upperElecBot )
    {
      retVal = 1;
    }
  else
    {
      retVal = 2;
    }
  return retVal;
}
