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
	void sample(Frame& a_frame);
	/// Normalizes the rdf by bin volume.
	void normalize();
	/// Prints the rdf to stdout.
	void print();
	/// Get number of pairs for which pair correlations must be computed.
	const unsigned int getNumPairs() const;
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
	/// Get address of first element in RDF
	double*** getRDFAddressLayers(int i);
	/// Get address of first element in 3d RDF
	double** getRDFAddressLayers(int i, int j);
	/// Get address of first element in 4d RDF
	double* getRDFAddressLayers(int i, int j, int k);
        /// Set value of element in RDF (for DEBUGGING ONLY)
        void setRDFLayerClosestValue(int a_layer, int a_bin, int a_pair, int a_closest, double a_setVal);

private:
	System m_system;
	vector<vector< double > > m_rdf;
	vector<vector<vector<double > > > m_rdfLayer;
	vector<vector<vector<vector<double > > > > m_rdfLayerClosest;
	vector<int > m_pairCounter;
	double m_maxDist;
	int m_numLayers;
	int m_numPairs;
	int m_numBins;
	double m_binSize;
	/// Puts pair distance into a bin.
	void binPairDistance(double a_distance, unsigned int a_pair);
	/// Puts pair distance into a bin based on the layer that the distance is contained in.
	void binPairDistance(double a_distance, unsigned int a_pair, unsigned int a_firstLayer, unsigned int a_secondLayer);
	/// Puts pair distance into a bin based on layer and reference species.
	void binPairDistance(double a_distance, unsigned int a_pair, unsigned int a_whichClosest, unsigned int a_firstLayer, unsigned int a_secondLayer);
	/// Increment count of pairs in standard RDF.
	void incrementCounter(unsigned int a_pair);
};

const char* RDFWrite(RDF* a_rdf, const char* a_filename);
const char* RDFWriteLayers(RDF* a_rdf, const char* a_filename);

#endif
