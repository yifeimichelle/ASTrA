#include "system.h"
#include "atom.h"
#include "frame.h"
#include "radialDistribution.h"
#include "writer.h"
#include <iostream>
#include <array>
using namespace std;

int main(int argc, char** argv)
{
    if(argc != 2)
    {
        cout << "this program takes one argument that is the input file ";
        cout << endl;
        return 1;

    }

    // read inputs
    string inputFile(argv[1]);
    System system(inputFile); // initialize system
    RDF rdf(system);
    Frame frame(system); // initialize trajectory frame reader
    for (int frameCounter = 0; frameCounter<system.getNumFrames(); frameCounter++)
      {
	frame.readStep();
	rdf.sample(frame);
	frame.clearFrame();
      }
    // normalize RDF
    rdf.normalize();

    // test indexing
    //double* rdfPtr = rdf.getRDFAddressLayers(0,0,0);
    //int idx;
    //cout << "testing RDF indexing ..." << endl;
    //rdf.setRDFLayerClosestValue(0,0,0,1,0.5);
    //cout << "0 0 0 1 value: " << rdf.getRDFLayerClosestElement(0,0,0,1) << endl;
    //idx = 1;
    //cout << "direct indexing " << rdfPtr[idx] << endl << endl;
    //rdf.setRDFLayerClosestValue(0,0,1,1,0.8);
    //cout << "0 0 1 1 value: " << rdf.getRDFLayerClosestElement(0,0,1,1) << endl;
    //idx = 1*2+1;
    //cout << "direct indexing " << rdfPtr[idx] << endl << endl;
    //rdf.setRDFLayerClosestValue(1,1,1,1,0.75);
    //cout << "1 1 1 1 value: " << rdf.getRDFLayerClosestElement(1,1,1,1) << endl;
    //idx = 1*(500*rdf.getNumPairs()*2) + 1*rdf.getNumPairs()*2 + 1*2 + 1;
    //cout << "direct indexing " << rdfPtr[idx] << endl << endl;
    //rdf.setRDFLayerClosestValue(2,10,5,1,1.5);
    //cout << "2 10 5 1 value: " << rdf.getRDFLayerClosestElement(2,10,5,1) << endl;
    //idx = 2*(500*rdf.getNumPairs()*2) + 10*rdf.getNumPairs()*2 + 5*2 + 1;
    //cout << "direct indexing " << rdfPtr[idx] << endl << endl;

    // print to a file
    RDFWrite(&rdf, "rdf");
    RDFWriteLayers(&rdf, "rdf");
}
