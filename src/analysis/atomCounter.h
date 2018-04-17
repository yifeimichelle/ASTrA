#ifndef _ATOMCOUNTER_H_
#define _ATOMCOUNTER_H_
#include "frame.h"
#include "system.h"
#include "atom.h"
#include <vector>
#define DENSITYCONV

using namespace std;
/// A class calculating atom counts and densities.
/**
 * This class stores number of atoms in each region (electrode/liquid), atom number and density profile across the system, and collective variables associated with atom counts.
 */

class AtomCounter
{
 public:
  AtomCounter();
  /// Constructor
  AtomCounter(System& a_system);
  /// Calculates COMs and computes density and number-of-atoms CVs given the current frame.
  void sample(const Frame& a_frame);
  /// Normalizes the density profiles and CVs.
  void normalize();
  /// Prints the density and number of atoms to stdout.
  void print();
  /// Returns the number of bins.
  const int getNumBins() const;
  /// Returns the bin size in length units.
  const int getBinSize() const;
  /// Bin atom.
  void binAtom(array<double, DIM>& a_position, int& a_molecType, int& a_molecMember);
  /// Bin molecule COM.
  void binElectrolyteCOM(array<double, DIM>& a_position, int& a_electrolyteID);
  void countElectrolyteInLayer(array<double, DIM>& a_position,  vector<array<int, 3> >& a_IonsInLayer, int& a_electrolyteID);

 private:
  System m_system;
  // Center of masses for this timestep.
  vector<array<double, DIM > > m_COMs;
  // Counter of all atoms in bins.
  vector<array<int, MAX_NUM_TYPES > > m_numAtomsProfile;
  // Counter of ion and solvent COMs in bins.
  vector<array<int, 3 > > m_numIonsProfile;
  // Counter of ion and solvent COMs in layers (anode, cathode, liquid).
  vector<array<int, 3 > > m_avgIonsInLayer;
  // Collective variables for atom and COM counts.
  vector<int > m_excessAnionsInCathode, m_excessCationsInAnode;
  vector<double > m_chargMechParam;
  int m_numBins;
  double m_binSize;
};

const char* ACWrite(AtomCounter* a_ac, const char* a_filename);
const char* ACWriteLayers(AtomCounter* a_ac, const char* a_filename);

#endif
