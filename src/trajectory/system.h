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
	/// Return number of frames
	const unsigned int getNumFrames() const;
	/// Return box dims
	const array<double, DIM > getBoxDims() const;
	/// Return whether the selected dimemsion is periodic
	const unsigned int isPeriodic(int i) const;
	/// Get number of layers (currently always returns 3, assuming anode, bulk, cathode layers)
	const unsigned int getNumLayers() const;
	/// Identify the layer that a set of coordinates is in
	const unsigned int getLayer(array<double, DIM>& a_position) const;
 private:
        string m_trajFile;
        unsigned int m_numFrames;
        unsigned int m_numMolecTypes;
	unsigned int m_numAtomTypes;
        unsigned int m_stepInterval;
	unsigned int m_numPairs;
	unsigned int m_cationID;
	unsigned int m_anionID;
	unsigned int m_anodeID;
	unsigned int m_cathodeID;
	unsigned int m_lowerElecID;
	unsigned int m_upperElecID;
	unsigned int m_anodeIsLower;
	unsigned int m_solventID;
	unsigned int m_capID;
	unsigned int m_boolWithSolvent;
	unsigned int m_boolWithCap;
	vector<vector<unsigned int> > m_typeAtomIndices;
	vector<pair<unsigned int, unsigned int> > m_rdfPairs;
        float m_stepTime;
        float m_frameTime;
	double m_lowerElecTop;
	double m_upperElecBot;
        array<double, DIM > m_boxDims;
        array<unsigned int, DIM > m_boxPeriodic;
        array<unsigned int, MAX_NUM_TYPES > m_numAtomsMolec;
        array<unsigned int, MAX_NUM_TYPES > m_numMolecs;
	unsigned int getAtomType(unsigned int a_molecType, unsigned int a_molecAtom);

};
#endif
