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
    //system.printTypeAtomIndices(); // for debugging
    int numAtoms = system.getNumAtoms();
    Frame frame(system); // initialize trajectory frame reader
    for (int frameCounter = 0; frameCounter<system.getNumFrames(); frameCounter++)
      {
	frame.readStep();
	array<double, DIM> position;
	position = frame.getAtom(672).getPosition();
	cout << "step " << frame.getStepNum() << ", position atom " << 672 << ": " << position[0] << " " << position[1] << " " << position[2] << endl;
	position = frame.getAtom(673).getPosition();
	cout << "step " << frame.getStepNum() << ", position atom " << 673 << ": " << position[0] << " " << position[1] << " " << position[2] << endl;
	//rdf.sample(frame);
	cout << frame.computeDistance(0,0) << endl;
	cout << frame.computeDistance(0,100) << endl;
	cout << frame.computeDistance(8000,1000) << endl;
	cout << frame.computeDistance(0,1) << endl;
	cout << frame.computeDistance(672,673) << endl;
	frame.clearFrame();
      }
    // normalize RDF
    // print to a file
}
