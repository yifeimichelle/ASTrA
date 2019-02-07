#ifndef _ATOMCOUNTER_H_
#define _ATOMCOUNTER_H_
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
    /// Samples zero-potential part of trajectory.
    void sampleZP(Frame& a_frame);
    /// Samples skipped part of trajectory, counting COMs and ions in layer vs. time only.
    void sampleSkip(Frame& a_frame);
    /// Normalizes the density profiles and CVs.
    void normalize();
    /// Normalizes the density profiles and CVs from the zero-potential, zero-Q run.
    void normalizeZP();
    /// Prints the number of atoms to stdout.
    void print();
    /// Prints density profile
    void printDensity();
    /// Returns the number of bins.
    const int getNumBins() const;
    /// Returns the bin size in length units.
    const double getBinSize() const;
    /// Returns the number of atom types in system.
    const int getNumAtomTypes();
    /// Returns the number of ion types (anion, cation, solvent) in the system.
    const int getNumIonTypes();
    /// Returns the number of layers in the system (typically 3: anode, cathode, liquid).
    const int getNumLayers() const;
    /// Gets the address of the first element in the vector of binned atom counts.
    double* getACAtomsAddress(int i);
    /// Gets the address of the first element in the vector of binned density.
    double* getACDensityAddress(int i);
    /// Gets the address of the first element in the vector of binned ion counts.
    double* getACIonsAddress(int i);
    /// Gets the address of the first element in the vector of ion counts per layer.
    double* getACIonsLayersAddress(int i);
    /// Gets the address of the element in the array of ions per layer per time.
    double* getACIonsLayersTimeAddress(int i, int j);
    /// Gets address of the first element in the charging mechanism CV.
    double* getChargingParamAddress(int i);
    /// Returns system.
    const System& getSystem() const;
    /// Returns interval for saving frames in time-sequence data.
    const int getSaveFrameInterval() const;
    /// Returns total number of saved frames (total frames)
    const int getNumSavedFrames() const;
    /// Returns number of saved frames for CV (only const-potential, const-charge frames)
    const int getNumCVFrames() const;
    /// Returns vector of COMs.
    //const vector<array<double, DIM > > getCOMs() const;

  private:
    /// Computes the charging mechanisms parameter from Forse 2016 "New perspectives on the charging mechanisms of supercapacitors".
    array<double, 2> computeChargingParam(vector<array<int, NUM_ION_TYPES> >& a_ionsInLayer);
    /// Bin and layer atom, and add to density profile, during zero-P, zero-Q run.
    void binZPAtom(Frame& a_frame, int& a_atomIndex, array<double, DIM>& a_position, int& a_molecType, int& a_molecMember, double& a_mass, int& a_isElectrolyte);
    /// Bin and layer electrolyte COM, during zero-P, zero-Q run.
    void binZPElectrolyteCOM(Frame& a_frame, int& a_molecIndex, array<double, DIM>& a_position, int& a_molecType, vector<array<int, 3> >& a_IonsInLayer,  int& a_electrolyteID, int& a_isElectrolyte);
    /// Compute number of ions in layers during skipped frames of run.
    void binSkipElectrolyteCOM(Frame& a_frame, int& a_molecIndex, array<double, DIM>& a_position, int& a_molecType, vector<array<int, 3> >& a_IonsInLayer,  int& a_electrolyteID, int& a_isElectrolyte);
    /// Bin and layer atom, and add to density profile.
    void binAtom(Frame& a_frame, int& a_atomIndex, array<double, DIM>& a_position, int& a_molecType, int& a_molecMember, double& a_mass, int& a_isElectrolyte);
    /// Bin and layer electrolyte COM.
    void binElectrolyteCOM(Frame& a_frame, int& a_molecIndex, array<double, DIM>& a_position, int& a_molecType, vector<array<int, 3> >& a_IonsInLayer,  int& a_electrolyteID, int& a_isElectrolyte);
    System m_system;
    /// Center of masses for this timestep.
    vector<array<double, DIM > > m_COMs;
    /// Collective variables for atom and COM counts, during entire trajectory (Z-P/Q and C-P/Q).
    vector<vector<array<double, NUM_ION_TYPES > > > m_numIonsInLayerTime;

    /// Counter of all atoms in bins.
    vector<array<double, MAX_NUM_TYPES > > m_numAtomsProfile;
    /// Stores density (from atom profile).
    vector<double > m_densityProfile;
    /// Counter of ion and solvent COMs in bins.
    vector<array<double, NUM_ION_TYPES > > m_numIonsProfile;
    /// Counter of ion and solvent COMs in layers (anode, cathode, liquid).
    vector<array<double, NUM_ION_TYPES > > m_avgIonsInLayer;

    /// Counter of all atoms in bins, during zero-P, zero-Q run.
    vector<array<double, MAX_NUM_TYPES > > m_numZPAtomsProfile;
    /// Stores density (from atom profile), during zero-P, zero-Q run.
    vector<double > m_ZPdensityProfile;
    /// Counter of ion and solvent COMs in bins, during zero-P, zero-Q run.
    vector<array<double, NUM_ION_TYPES > > m_numZPIonsProfile;
    /// Counter of ion and solvent COMs in layers (anode, cathode, liquid), during zero-P, zero-Q run.
    vector<array<double, NUM_ION_TYPES > > m_avgZPIonsInLayer;
    vector<int > m_excessAnionsInCathode, m_excessCationsInAnode;
    vector<array<double, 2> > m_chargingParam;
    int m_saveFrameEvery;
    int m_numSavedFrames;
    int m_numBins;
    double m_binSize;
    int m_numAtomTypes;
    int m_numLayers;
    double m_zLo;
};

const char* ACWriteAtomCounts(AtomCounter* a_ac, const char* a_filename);
const char* ACWriteDensity(AtomCounter* a_ac, const char* a_filename);
const char* ACWriteIons(AtomCounter* a_ac, const char* a_filename);
const char* ACWriteIonsInLayers(AtomCounter* a_ac, const char* a_filename);
const char* ACWriteCollectiveVars(AtomCounter* a_ac, const char* a_filename);
const char* ACWriteIonsInLayersTime(AtomCounter* a_ac, const char* a_filename);


#endif
