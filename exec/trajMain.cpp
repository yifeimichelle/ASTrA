#include "system.h"
#include "atom.h"
#include "frame.h"
#include "radialDistribution.h"
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
    int numAtoms = system.getNumAtoms();
    Frame frame(system); // initialize trajectory frame reader
    for (int frameCounter = 0; frameCounter<system.getNumFrames(); frameCounter++)
      {
	frame.readStep();
	rdf.sample(frame);
	frame.clearFrame();
      }
    // normalize RDF
    rdf.normalize();
    //rdf.print();
    // print to a file
    RDFWrite(&rdf, "rdf");
    RDFWriteLayers(&rdf, "rdf");

}
