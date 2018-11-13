#ifndef _FRAME_H_
#define _FRAME_H_

#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include "system.h"
#include "atom.h"
#include <list>

using namespace std;
/// A class describing a single frame (snapshot) of a trajectory
/**
 * This class stores atomic positions
 */

class Frame
{
    public:
        Frame();
        /// Constructor for initial frame.
        Frame(System& a_system);
	/// Reads a constant-P or -Q step. Stores atom positions and increments m_stepNum.
        void readStep();
	/// Reads a constant-P or -Q step. Stores atom positions and increments m_stepNum.
        void readStep(int a_every);
	/// Skips a step. Doesn't store anything. Only increments m_totalStepNum.
	void skipStep();
	/// Skips a step. Doesn't store anything. Only increments m_totalStepNum.
	void skipStep(int a_every);
	/// Reads a zero-potential, zero-charge step. Stores atom positions and increments m_zpStepNum.
	void readZPStep();
	/// Reads a zero-potential, zero-charge step. Stores atom positions and increments m_zpStepNum.
	void readZPStep(int a_every);
	/// Reads fluctuating charges for electrodes.
	void readCharges();
	/// Skips fluctuating charges for electrodes.
	void skipCharges();
	/// Returns total step number (counting from beginning of trajectory and not skipping any steps).
	const unsigned int getTotalStepNum() const;
	/// Returns step number of the current frame in the constant-P or -Q window.
	const unsigned int getStepNum() const;
	/// Returns step number of the current frame in the zero-P or -Q window.
	const unsigned int getZPStepNum() const;
	/// Sets all atomic coordinates to -1.
        void clearFrame();
	/// Returns a reference to a single atom in a frame.
	const Atom& getAtom(int a_atomIndex) const;
	/// Returns a reference to a single atom in a frame.
	Atom& getAtom(int a_atomIndex);
	/// Returns a reference to a single molecule in a frame.
	const Atom& getMolec(int a_molecIndex) const;
	/// Returns a reference to a single molecule in a frame.
	Atom& getMolec(int a_molecIndex);
	/// Returns a reference to a single atom in a frame.
	const Atom& getAtomOfMolec(int a_molecIndex) const;
	/// Computes the distance between two atoms.
	const double computeDistance(int a_i, int a_j) const;
	/// Computes the distance between two molecules.
	const double computeMolecDistance(int a_i, int a_j) const;
	/// Identify the layer an atom is in
	const unsigned int getLayerOf(unsigned int a_index) const;
	/// Identify the layer a molecule is in
	const unsigned int getLayerOfMolec(unsigned int a_index) const;
	/// Get COMs of molecules from AtomCounter
	void setCOMs(vector<array<double, DIM > > a_COMs);
	/// Assign atom to layer, during zero-P, zero-Q run.
	void assignZPAtomToLayer(unsigned int a_index, unsigned int a_type, unsigned int a_layer);
	/// Assign electrolyte COM to layer, during zero-P, zero-Q run.
	void assignZPIonToLayer(unsigned int a_index, unsigned int a_type, unsigned int a_layer);
	/// Assign atom to layer.
	void assignAtomToLayer(unsigned int a_index, unsigned int a_type, unsigned int a_layer);
	/// Assign electrolyte COM to layer.
	void assignIonToLayer(unsigned int a_index, unsigned int a_type, unsigned int a_layer);
	/// Print atoms in layer.
	void printAtomsInLayer(unsigned int a_layer);
	/// Print atoms in layer other way (for debugging only)
	void printAtomsInLayerCheck(unsigned int a_layer);
	/// Print molecules in layer.
	void printMolecsInLayer(unsigned int a_layer);
	/// Get pointer to atoms in layer.
	vector<int>* getAtomsInLayer(int a_layerIdx) const;
	/// Get pointer to molecule COMs in layer.
	vector<int>* getMoleculesInLayer(int a_layerIdx) const;
	/// Computes the degree of confinement for the given frame
	void computeDegreeOfConfinement();
	/// Gets number of atoms in layer.
	const unsigned int getCurrentNumAtomsInLayer(int a_layerIdx, int a_molID) const;
	/// Gets number of molecules in layer.
	const unsigned int getCurrentNumMolecsInLayer(int a_layerIdx, int a_molID) const;
	/// Sums charges of selected molecule.
	const double sumCharges(const System& a_system, int a_molID) const;
     private:
	/// Set charges to charges read from input file.
	void setCharges(const System& a_system);
	System m_system;
	unsigned int m_every;
	unsigned int m_stepNum;
	unsigned int m_totalStepNum;
	unsigned int m_zpStepNum;
        unsigned int m_numAtoms;
        unsigned int m_numMolecules;
	unsigned int m_numFluctuatingCharges;
        vector<Atom > m_atoms;
	vector<Atom > m_COMs;
	// record of which layer atoms and molecule COMs are in
	// access by ARRAY[layer][type] -> list it
	array<array<vector<int >, MAX_NUM_TYPES> , NUM_LAYERS > m_atomLayers;
	array<array<vector<int >, MAX_NUM_TYPES> , NUM_LAYERS > m_COMLayers;
	array<array<vector<int >, MAX_NUM_TYPES> , NUM_LAYERS > m_ZPatomLayers;
	array<array<vector<int >, MAX_NUM_TYPES> , NUM_LAYERS > m_ZPCOMLayers;
        ifstream m_traj;
	ifstream m_chg;
};
#endif
