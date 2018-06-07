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
    cout << "number of AC bins: " << ac.getNumBins() << endl;
    Frame frame(system); // initialize trajectory frame reader

    cout << "Reading trajectory ..." << endl;
    cout << "   box x and y: " << system.getBoxDim(0) << " " << system.getBoxDim(1) << endl;

    // Skip frames
    if (system.getNumZPFrames() > 0)
      {
	for (int i=0; i<system.getNumSkipFrames(); i++)
	  {
	    frame.skipStep();
	    ac.sampleSkip(frame);
	    if (READ_CHARGE_FILE)
	      {
		frame.skipCharges(); //!!
	      }
	  }

	// Read zero-potential, zero-charge frames
	cout << "Analyzing zero-P, zero-Q run of " << system.getNumZPFrames() << " steps..." << endl;
	for (unsigned int frameCounter = 0; frameCounter<system.getNumZPFrames(); frameCounter++)
	  {
	    frame.readZPStep();
	    if (READ_CHARGE_FILE)
	      {
		frame.skipCharges(); //!!
	      }

	    if (frame.getZPStepNum() % int(ceil(system.getNumTotalFrames()/10.0)) == 0)
	      {
		cout << frame.getZPStepNum() << endl;
	      }

	    ac.sampleZP(frame);
	    frame.clearFrame();
	  }
	ac.normalizeZP();
      }

    // Skip frames (potential or charge turned on)
    for (int i=0; i<system.getNumSkipFrames(); i++)
      {
	frame.skipStep();
	ac.sampleSkip(frame);
      }

    // Read constant-potential or constant-charge frames
    cout << "Analyzing constant-P or -Q run of " << system.getNumFrames() << " steps..." << endl;
    for (unsigned int frameCounter = 0; frameCounter<system.getNumFrames(); frameCounter++)
      {
	// read in step of trajectory
	frame.readStep();
	if (READ_CHARGE_FILE)
	  {
	    frame.readCharges(); //!!
	  }

	if ( frame.getStepNum() % int(ceil(system.getNumTotalFrames()/10.0)) == 0)
	  {
	    cout << frame.getStepNum() << endl;
	  }

	// sample routines
	ac.sample(frame);

	// clear frame memory
	frame.clearFrame();
      }
    // normalize atom count / density profile
    ac.normalize();

    // print to a file
    ACWriteAtomCounts(&ac, "atoms");
    ACWriteDensity(&ac, "density");
    ACWriteIons(&ac, "ions");
    ACWriteIonsInLayers(&ac, "layers");
    ACWriteIonsInLayersTime(&ac, "numionslayers");
    ACWriteCollectiveVars(&ac, "ionCV");
}
