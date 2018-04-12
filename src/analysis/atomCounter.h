#include "frame.h"
#include "system.h"
#include "atom.h"
#include <vector>

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
  void sample(Frame& a_frame);
  /// Normalizes the density profiles and CVs.
  void normalize();
  /// Prints the density and number of atoms to stdout.
  void print();
  /// Returns the number of bins.
  const int getNumBins() const;
  /// Returns the bin size in length units.
  const int getBinSize() const;
  /// Bin atom.
  void binAtom(array<double, DIM> a_position, int a_molecType, int a_molecMember);
  /// Bin molecule COM.
  void binElectrolyteCOM(array<double, DIM> a_position, int a_electrolyteID);

 private:
  System m_system;
  // Center of masses for this timestep.
  vector<array<double, DIM > > m_COMs;
  // Counter of all atoms in cell.
  vector<array<int, MAX_NUM_TYPES > > m_numAtomsProfile;
  // Counter of ion and solvent COMs (max 3 for organic electrolyte).
  vector<array<int, 3 > > m_numIonsProfile;
  int m_numBins;
  double m_binSize;
  /// Compute COMs
  void computeCOMs(Frame& a_frame);
  /// Compute density profile
  void computeDensity(Frame& a_frame);
};
