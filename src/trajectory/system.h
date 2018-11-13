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
#define NUM_ION_TYPES 3
#define NUM_LAYERS 3
#ifndef READ_CHARGE_FILE
#define READ_CHARGE_FILE 0
#endif

using namespace std;
/// A class describing the system (either for MD or MC).
/**
 * This class stores atom types and indices, box and periodic boundary specifications.
 */

class System
{
    public:
        /// Default constructor.
        System();
        /// Constructor from input file.
        System(const string& a_inputFile);
	/// New constructor.
	//System(const string& a_inputFile, int a_
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
	/// Get total number of molecules in the system (including mobile and immobile)
	const unsigned int getNumMolecules() const;
	/// Return real time interval between frames;
        const double getFrameTime() const;
	/// Returns the name of the trajectory file.
        const string& getTrajFile() const;
	/// Returns the name of the fluctuating electrode charges file.
        const string& getChargesFile() const;
	/// Returns number of atoms of a particular atom type (each atom in a molecule counts separately).
	const unsigned int getNumOfType(unsigned int a_type) const;
	/// Returns index of individual atom, given atom type and index into all atoms of that atom type.
	const unsigned int getIndexOfType(unsigned int a_type, unsigned int a_idx) const;
	/// Get index of molecule for COMs.
	const unsigned int getMolecIndex(unsigned int a_type, unsigned int a_member) const;
	/// Return one-dimensional atom type index
	const unsigned int getAtomType(unsigned int a_molecType, unsigned int a_molecAtom) const;
	/// Returns number of RDF atom-atom pairs to calculate, supplied in input.
	const unsigned int getNumPairs() const;
	/// Returns pair of atom types from list of RDF atom-atom pairs.
	const pair<unsigned int, unsigned int >& getPairCorrelation(unsigned int a_pair) const;
	/// Returns number of RDF molecule-molecule pairs to calculate, supplied in input.
	const unsigned int getNumMolecPairs() const;
	/// Returns pair of molecule types from list of RDF molecule-molecule pairs.
	const pair<unsigned int, unsigned int > getMolecPairCorrelation(unsigned int a_pair) const;
	/// Returns coordination number cutoff for RDF molecule-molecule pairs.
	const double getMolecPairCutoff(unsigned int a_pair) const;
	/// Return box dims.
	const array<double, DIM >& getBoxDims() const;
	/// Return box size in specified dimension.
	const double getBoxDim(int a_dim) const;
        /// Return whether the box is symmetric over the z=0 axis.
        const unsigned int isZSymmetrized() const;
	/// Return lower z boundary
	const double getZLo() const;
	/// Return whether the selected dimemsion is periodic.
	const unsigned int isPeriodic(int i) const;
	/// Get number of layers (currently always returns 3, assuming anode, bulk, cathode layers).
	const unsigned int getNumLayers() const;
	/// Return number of frames in the trajectory.
	const unsigned int getNumFrames() const;
	/// Return total number of frames in the trajectory.
	const unsigned int getNumTotalFrames() const;
	/// Return number of zero-potential, uncharged frames in the trajectory.
	const unsigned int getNumZPFrames() const;
	/// Return number of skipped equilibration frames in the trajectory.
	const unsigned int getNumSkipFrames() const;
	/// Get number of species in electrolyte (2 if IL, 3 if organic electrolyte)
	const unsigned int getNumElectrolyteSpecies() const;
	/// Get total number of molecules of electrolytes
	const unsigned int getNumElectrolyteMolecs() const;
	/// Identify the layer that a set of coordinates is in.
	const unsigned int getLayer(array<double, DIM>& a_position) const;
	/// Return masses of atoms in specified molecule type.
	const array<double, MAX_MEMBERS_PER_MOLEC > getMassesOfType(int a_type) const;
	/// Return charges of atoms in specified molecule type.
	const array<double, MAX_MEMBERS_PER_MOLEC > getChargesOfType(int a_type) const;
	/// Returns whether molecule is an electrolyte component and stores ID of component (cation, anion, solvent) in second argument.
	unsigned int isElectrolyte(int a_molecType, int* a_electrolyteID) const;
	/// Inserts layers into array pointer
	void getLayerUpperBounds(int a_numLayers, double* a_layers) const;
	/// Returns whether anode is the "lower" electrode in the system.
	unsigned int isAnodeLower() const;
	/// Returns whether ID is cathode molecule.
	unsigned int isCathode(unsigned int a_molID) const;
	/// Returns whether ID anode molecule.
	unsigned int isAnode(unsigned int a_molID) const;
	/// Returns whether or not charge file is to be read.
	unsigned int hasChargeFile() const;
	/// Returns anode ID.
	unsigned int getAnodeID() const;
	/// Returns cathode ID.
	unsigned int getCathodeID() const;
	/// Returns whether molecule is anion.
	unsigned int isAnion(unsigned int a_molID) const;
	/// Returns whether molecule is cation.
	unsigned int isCation(unsigned int a_molID) const;
	/// Returns first atom index of molecule.
	const unsigned int getFirstAtomOfMolec(unsigned int a_molecIndex) const;
	/// Returns read interval for trajectory steps
	const unsigned int getReadFrameEvery() const;
 private:
	void readInput(const string& a_inputFile);
	void setInput();
	template <typename T>
	  void getInput(T* a_value, int a_col);
	template <typename T1, typename T2>
	  void getInputs2(T1* a_value1, T2* a_value2);
	template <typename T1, typename T2, typename T3>
	  void getInputs3(T1* a_value1, T2* a_value2, T3* a_value3);
	template <typename T1, typename T2, typename T3, typename T4>
	  void getInputs4(T1* a_value1, T2* a_value2, T3* a_value3, T4* a_value4);
	template <typename T1, typename T2>
	  void getInputs2(T1* a_value1, T2* a_value2, int a_offset);
	template <typename T1, typename T2, typename T3>
	  void getInputs3(T1* a_value1, T2* a_value2, T3* a_value3, int a_offset);
	vector<vector<string > > m_inputs;
	void nextRow();
	unsigned int m_inputRow;
	string m_trajFile;
	string m_chgFile;
	unsigned int m_totalFrames;
        unsigned int m_numFramesInclSkip; // Number of constant-P or constant-Q "production" frames.
        unsigned int m_numFrames; // Number of constant-P or constant-Q "production" frames.
        unsigned int m_zpFramesInclSkip;
        unsigned int m_zpFrames;
        unsigned int m_skipFrames;
        unsigned int m_numMolecTypes;
	unsigned int m_numAtoms;
	unsigned int m_numAtomTypes;
        unsigned int m_stepInterval;
	unsigned int m_numPairs;
	unsigned int m_numMolecPairs;
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
	unsigned int m_numElectrolyteMolecs;
	unsigned int m_numMolecules;
	unsigned int m_readFrameEvery;
        unsigned int m_zSymmetrized;
	double m_zLo;
	vector<int > m_firstAtomOfMolec;
	array<array<double , MAX_MEMBERS_PER_MOLEC >, MAX_NUM_TYPES > m_masses;
	array<array<double , MAX_MEMBERS_PER_MOLEC >, MAX_NUM_TYPES > m_charges;
	vector<vector<unsigned int > > m_typeAtomIndices;
	vector<pair<unsigned int, unsigned int > > m_molecMembersOfType;
	vector<pair<unsigned int, unsigned int > > m_rdfPairs;
	vector<pair<unsigned int, unsigned int > > m_rdfMolecPairs;
	vector<double > m_rdfMolecCutoffs;
        float m_stepTime;
        float m_frameTime;
	unsigned int m_readFluctuatingCharge;
	double m_lowerElecTop;
	double m_upperElecBot;
        array<double, DIM > m_boxDims;
        array<unsigned int, DIM > m_boxPeriodic;
        array<unsigned int, MAX_NUM_TYPES > m_numMembersMolec;
        array<unsigned int, MAX_NUM_TYPES > m_numMolecs;
	vector<string > lineToString(char* a_inputline, string& a_delimiter);
	vector<string > readNextLine(char* a_inputline, string& a_delimiter);
};
#endif
