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
	/// Reads a step. Stores atom positions and increments m_stepNum.
        void readStep();
	/// Returns step number of the current frame.
	const unsigned int getStepNum() const;
	/// Sets all atomic coordinates to -1.
        void clearFrame();
	/// Returns a reference to a single atom in a frame.
	const Atom& getAtom(int a_atomIndex) const;
	/// Returns a reference to a single atom in a frame.
	Atom& getAtom(int a_atomIndex);
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
	/// Assign atom to layer.
	void assignAtomToLayer(unsigned int a_index, unsigned int a_type, unsigned int a_layer);
	/// Assign electrolyte COM to layer.
	void assignIonToLayer(unsigned int a_index, unsigned int a_type, unsigned int a_layer);
	/// Print atoms in layer.
	void printAtomsInLayer(unsigned int a_layer);
	/// Print atoms in layer other way (for debugging only)
	void printAtomsInLayerCheck(unsigned int a_layer);
	/// Get pointer to atoms in layer.
	list<int>* getAtomsInLayer(int a_layerIdx) const;
     private:
	System m_system;
        unsigned int m_stepNum;
        unsigned int m_numAtoms;
        unsigned int m_numMolecs;
        vector<Atom > m_atoms;
	vector<Atom > m_COMs;
	// record of which layer atoms and molecule COMs are in
	// access by ARRAY[layer][type] -> list it
	array<array<list<int >, MAX_NUM_TYPES> , NUM_LAYERS > m_atomLayers;
	array<array<list<int >, MAX_NUM_TYPES> , NUM_LAYERS > m_COMLayers;
        ifstream m_traj;
};
#endif
