#ifndef _RADIALDISTRIBUTION_H_
#define _RADIALDISTRIBUTION_H_
#include "frame.h"
#include "system.h"
#include "atom.h"
#include <vector>

using namespace std;
/// A class to contain radial distribution (pair correlation) functions and associated calculations.
/**
 * This class stores parameters and data for the radial distribution functions. RDFs calculated: Standard RDF, LRDF (layer-specific RDF), and CLRDF (closest-neighbor layer-specific RDF).
 */

class RDF
{
    public:
        RDF();
        /// Constructor 
        RDF(System& a_system);
	/// Computes radial distribution functions given the current frame.
	void sample(const Frame& a_frame);
	/// Computes atom-atom radial distribution functions given the current frame.
	void sampleAtoms(const Frame& a_frame);
	/// Computes molecule-molecule COM radial distribution functions given the current frame.
	void sampleMolecules(const Frame& a_frame);
	/// OLD: Computes radial distribution functions given the current frame.
	void sampleOld(Frame& a_frame);
	/// Samples zero-potential part of trajectory.
	void sampleZP(Frame& a_frame);
	/// Normalizes the rdf by bin volume.
	void normalize();
	/// Prints the rdf to stdout.
	void print();
	/// Get number of pairs for which pair correlations must be computed.
	const unsigned int getNumPairs() const;
	/// Get number of molecule-molecule COM pairs for which pair correlations must be computed.
	const unsigned int getNumMolecPairs() const;
	/// Get number of bins over which pair correlations are computed.
	const unsigned int getNumBins() const;
	/// Get size of bins (in length units).
	const double getBinSize() const;
	/// Get value of element in standard RDF.
	const double getRDFElement(int a_distance, int a_pair) const;
	/// Get value of element in LRDF.
	const double getRDFLayerElement(unsigned int a_layer, int a_distance, int a_pair) const;
	/// Get value of element in CLRDF.
    	const double getRDFLayerClosestElement(unsigned int a_layer, int a_distance, int a_pair, int a_whichClosest) const;
	/// Get number of layers in system (is 3 by default with anode, cathode, liquid).
	const unsigned int getNumLayers() const;
	/// Get address of first element in RDF
	double* getRDFAddress(int i);
	/// Get address of first element in 3d RDF
	double* getRDFAddressLayers(int i, int j);
	/// Get address of first element in 4d RDF
	double* getRDFAddressLayersClosest(int i, int j, int k);
	/// Get address of first element in molecule-molecule RDF
	double* getMolecRDFAddress(int i);
	/// Get address of first element in 3d molecule-molecule RDF
	double* getMolecRDFAddressLayers(int i, int j);
	/// Get address of first element in 4d molecule-molecule RDF
	double* getMolecRDFAddressLayersClosest(int i, int j, int k);
	/// Get address of first element in DoC
	double* getDoCAddress(int i);
	/// Set value of element in RDF (for DEBUGGING ONLY)
        void setRDFLayerClosestValue(int a_layer, int a_bin, int a_pair, int a_closest, double a_setVal);
	/// Clear data for the current frame
	void clearFrame();
	/// Compute degree of confinement
	void computeDegreeOfConfinement(const Frame& a_frame);

private:
	System m_system;
	vector<vector< double > > m_rdf;
	vector<vector<vector<double > > > m_rdfLayer;
	vector<vector<vector<vector<double > > > > m_rdfLayerClosest;
	vector<int > m_pairCounter;
	vector<vector< double > > m_rdfMolec;
	vector<vector<vector<double > > > m_rdfMolecLayer;
	vector<vector<vector<double > > > m_currentRDFMolecLayer;
	vector<vector<vector<vector<double > > > > m_rdfMolecLayerClosest;
	vector<int > m_pairMolecCounter;
	double m_maxDist;
	int m_numLayers;
	int m_numPairs;
	int m_numMolecPairs;
	int m_numBins;
	double m_binSize;
	// Constants for calculating DoC
	double m_Rj;
	double m_phi;
	double m_RcutDoC;
	vector<double > m_solidAngleFactor;
	vector<vector<double > > m_DoC;
	// end constants for calculating DoC
	/// Puts atom-atom pair distance into a bin based on layer.
	void binPairDistanceLayer(double a_distance, unsigned int a_pair, unsigned int a_layer);
	/// Puts atom-atom pair distance into a bin for the specified layer.
	void binPairDistanceClosestLayer(double a_distance, unsigned int a_pair, unsigned int a_whichClosest, unsigned int a_layer);
	/// Puts molecule-molecule COM pair distance into a bin based on the layer.
	void binMolecPairDistanceLayer(double a_distance, unsigned int a_pair, unsigned int a_layer);
	/// Puts molecule-molecule COM pair distance into a bin for the specified layer.
	void binMolecPairDistanceClosestLayer(double a_distance, unsigned int a_pair, unsigned int a_whichClosest, unsigned int a_layer);
	/// Increment count of atom-atom pairs in standard RDF.
	void incrementCounter(unsigned int a_pair);
	/// Increment count of molecule-molecule pairs in standard RDF.
	void incrementMolecCounter(unsigned int a_pair);
};

const char* RDFWrite(RDF* a_rdf, const char* a_filename);
const char* RDFWriteLayers(RDF* a_rdf, const char* a_filename);
const char* RDFWriteLayersClosest(RDF* a_rdf, const char* a_filename);
const char* RDFMolecWrite(RDF* a_rdf, const char* a_filename);
const char* RDFMolecWriteLayers(RDF* a_rdf, const char* a_filename);
const char* RDFMolecWriteLayersClosest(RDF* a_rdf, const char* a_filename);
const char* DoCWrite(RDF* a_rdf, const char* a_filename);

#endif
