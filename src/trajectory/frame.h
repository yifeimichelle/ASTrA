#ifndef _FRAME_H_
#define _FRAME_H_

#include <vector>
#include <string>
#include <cstring>
#include <fstream>
#include "system.h"
#include "atom.h"

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
	double computeDistance(int a_i, int a_j);
	/// Identify the layer an atom is in
	unsigned int getLayerOf(unsigned int a_index);
    private:
	System m_system;
        unsigned int m_stepNum;
        unsigned int m_numAtoms;
        vector<Atom > m_atoms;
        ifstream m_traj;
};
#endif
