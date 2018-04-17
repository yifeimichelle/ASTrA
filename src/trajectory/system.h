#ifndef _SYSTEM_H_
#define _SYSTEM_H_

#include <cstdio>
#include <cmath>
#include <cstring>
#include <array>
#include <vector>
#include <string>
#define MAX_NUM_TYPES 10
#define MAX_MEMBERS_PER_MOLEC 3

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
	/// Prints pairs specified in input to stdout (useful for debugging).
	void printPairCorrelations() const;
	/// Prints atom types and indices to stdout (useful for debugging).
	void printTypeAtomIndices() const;
	/// Returns total number of atoms in the system.
        const int getNumAtoms() const;
	/// Returns number of unique atom types.
	const int getNumAtomTypes() const;
	/// Returns number of molecule types in the system (e.g., ions, solvent, anode, cathode, cap).
	const int getNumMolecTypes() const;
	/// Returns number of molecules of specified type.
	const int getNumMolecsOfType(unsigned int a_molecType) const;
	/// Returns number of member atoms comprising a molecule.
	const int getNumMembersMolec(unsigned int a_molecType) const;
	/// Return real time of current frame;
        const float getFrameTime() const;
	/// Returns the name of the trajectory file.
        const string& getTrajFile() const;
	/// Returns number of atoms of a particular atom type (each atom in a molecule counts separately).
	const unsigned int getNumOfType(unsigned int a_type) const;
	/// Returns index of individual atom, given atom type and index into all atoms of that atom type.
	const unsigned int getIndexOfType(unsigned int a_type, unsigned int a_idx) const;
	/// Get index of molecule for COMs.
	const unsigned int getMolecIndex(unsigned int a_type, unsigned int a_member) const;
	/// Return one-dimensional atom type index
	const unsigned int getAtomType(unsigned int a_molecType, unsigned int a_molecAtom) const;
	/// Returns number of RDF pairs to calculate, supplied in input.
	const unsigned int getNumPairs() const;
	/// Returns pair of atom types from list of RDF pairs.
	const pair<unsigned int, unsigned int > getPairCorrelation(unsigned int a_pair) const;
	/// Return box dims.
	const array<double, DIM >& getBoxDims() const;
	/// Return box size in specified dimension.
	const double getBoxDim(int a_dim) const;
	/// Return whether the selected dimemsion is periodic.
	const unsigned int isPeriodic(int i) const;
	/// Get number of layers (currently always returns 3, assuming anode, bulk, cathode layers).
	const unsigned int getNumLayers() const;
	/// Return number of frames in the trajectory.
	const unsigned int getNumFrames() const;
	/// Get number of species in electrolyte (2 if IL, 3 if organic electrolyte)
	const unsigned int getNumElectrolyteSpecies() const;
	/// Identify the layer that a set of coordinates is in.
	const unsigned int getLayer(array<double, DIM>& a_position) const;
	/// Return masses of atoms in specified molecule type.
	const array<double, MAX_MEMBERS_PER_MOLEC > getMassesOfType(int a_type) const;
	/// Return charges of atoms in specified molecule type.
	const array<double, MAX_MEMBERS_PER_MOLEC > getChargesOfType(int a_type) const;
	/// Returns whether molecule is an electrolyte component and stores ID of component (cation, anion, solvent) in second argument.
	unsigned int isElectrolyte(int a_molecType, int* a_electrolyteID) const;
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
	unsigned int m_numElectrolyteSpecies;
	array<array<double , MAX_MEMBERS_PER_MOLEC >, MAX_NUM_TYPES > m_masses;
	array<array<double , MAX_MEMBERS_PER_MOLEC >, MAX_NUM_TYPES > m_charges;
	vector<vector<unsigned int > > m_typeAtomIndices;
	vector<pair<unsigned int, unsigned int > > m_molecMembersOfType;
	vector<pair<unsigned int, unsigned int > > m_rdfPairs;
        float m_stepTime;
        float m_frameTime;
	double m_lowerElecTop;
	double m_upperElecBot;
        array<double, DIM > m_boxDims;
        array<unsigned int, DIM > m_boxPeriodic;
        array<unsigned int, MAX_NUM_TYPES > m_numMembersMolec;
        array<unsigned int, MAX_NUM_TYPES > m_numMolecs;

};
#endif
