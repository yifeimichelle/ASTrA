#ifndef _SYSTEM_H_
#define _SYSTEM_H_

#include <cstdio>
#include <cmath>
#include <cstring>
#include <array>
#include <vector>
#include <string>
#define MAX_NUM_TYPES 10

using namespace std;
/// A class describing the system (either for MD or MC).
/**
 * This class stores atom types and indices, box and periodic boundary specifications.
 */

class System
{
    public:
        System();
        /// Constructor from input file.
        System(const std::string& a_inputFile);
        const int getNumAtoms() const;
        const float getFrameTime() const;
        const string& getTrajFile() const;
	void printPairCorrelations() const; // change this to print to a file
	void printTypeAtomIndices() const; // change this to print to a file
	const unsigned int getNumOfType(unsigned int a_type) const;
	const unsigned int getIndexOfType(unsigned int a_type, unsigned int a_idx) const;
	const unsigned int getNumPairs() const;
	const pair<unsigned int, unsigned int > getPairCorrelation(unsigned int a_pair) const;
	const unsigned int getNumFrames() const;
	const array<double, DIM > getBoxDims() const;
	const unsigned int isPeriodic(int i) const;
 private:
        string m_trajFile;
        unsigned int m_numFrames;
        unsigned int m_numMolecTypes;
	unsigned int m_numAtomTypes;
        unsigned int m_stepInterval;
	unsigned int m_numPairs;
	vector<vector<unsigned int> > m_typeAtomIndices;
	vector<pair<unsigned int, unsigned int> > m_rdfPairs;
        float m_stepTime;
        float m_frameTime;
        array<double, DIM > m_boxDims;
        array<unsigned int, DIM > m_boxPeriodic;
        array<unsigned int, MAX_NUM_TYPES > m_numAtomsMolec;
        array<unsigned int, MAX_NUM_TYPES > m_numMolecs;
	unsigned int getAtomType(unsigned int a_molecType, unsigned int a_molecAtom);

};
#endif
