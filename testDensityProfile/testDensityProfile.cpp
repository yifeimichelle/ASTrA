#include "system.h"
#include "atom.h"
#include "frame.h"
#include "atomCounter.h"
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
    AtomCounter ac(system);
    cout << "number of bins: " << ac.getNumBins() << endl;
    Frame frame(system); // initialize trajectory frame reader
    for (unsigned int frameCounter = 0; frameCounter<system.getNumFrames(); frameCounter++)
      {
	frame.readStep();
	array<double, DIM> testPosition = frame.getAtom(0).getPosition();
	cout << "read  : " << testPosition[0] << " " << testPosition[1] << " " << testPosition[2] << endl;
	int tmp=0;
	ac.binElectrolyteCOM(testPosition, tmp);
	ac.binAtom(testPosition, tmp, tmp);
	ac.sample(frame);
	testPosition = frame.getAtom(0).getPosition();
	cout << "sample: " << testPosition[0] << " " << testPosition[1] << " " << testPosition[2] << endl;
	ac.print();
	frame.clearFrame();
	testPosition = frame.getAtom(0).getPosition();
	cout << "clear : " << testPosition[0] << " " << testPosition[1] << " " << testPosition[2] << endl;
	cout << endl;
      }
    // normalize atom count / density profile
    ac.normalize();
    // print to a file
    ACWriteDensity(&ac, "density");
    ACWriteIons(&ac, "ions");
    ACWriteIonsInLayers(&ac, "layers");
}
