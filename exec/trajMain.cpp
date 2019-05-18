#include "system.h"
#include "atom.h"
#include "frame.h"
#include "radialDistribution.h"
#include "atomCounter.h"
#include <iostream>
#include <array>
using namespace std;

int main(int argc, char** argv)
{
  if(argc != 2) //!! Make it able to take additional arguments
  {
    cout << "this program takes one argument that is the input file ";
    cout << endl;
    return 1;
  }

  // read inputs
  string inputFile(argv[1]);
  System system(inputFile); // initialize system
  Frame frame(system); // initialize trajectory frame reader
  AtomCounter ac(system); // initialize atomcounter
  RDF rdf(system, ac); // initialize rdf with system AND atomcounter

  cout << "Reading trajectory ..." << endl;

  // Skip frames
  if (system.getNumZPFrames() > 0)
  {
    for (int i=0; i<system.getNumSkipFrames(); i++)
    {
      frame.skipStep();
      if (READ_CHARGE_FILE)
      {
        frame.skipCharges();
      }
      ac.sampleSkip(frame);
    }

    // Read zero-potential, zero-charge frames
    cout << "Analyzing zero-P, zero-Q run of " << system.getNumZPFrames() << " steps..." << endl;
    for (unsigned int frameCounter = 0; frameCounter<system.getNumZPFrames(); frameCounter++)
    {
      frame.readZPStep();
      if (READ_CHARGE_FILE)
      {
        frame.skipCharges();
      }

      if (frame.getZPStepNum() % int(ceil(system.getNumTotalFrames()/10.0)) == 0)
      {
        cout << frame.getZPStepNum() << endl;
        int timestep = frame.getTimestep();
        cout << "Current timestep is " << timestep << endl;
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
    if (READ_CHARGE_FILE)
    {
      frame.skipCharges();
    }
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
      int timestep = frame.getTimestep();
      cout << "Current timestep is " << timestep << endl;
    }

    // sample routines
    ac.sample(frame);
    rdf.sample(frame);
    rdf.computeDegreeOfConfinement(frame);

    // clear frame memory
    frame.clearFrame();
    rdf.clearFrame();
  }
  // normalize RDF
  ac.normalize();
  rdf.normalize(&ac);

  // print to stdout
  //ac.printDensity();
  //rdf.print();
  //ac.print();

  // print to a file
  RDFWrite(&rdf, "rdf");
  RDFWriteLayers(&rdf, "rdf");
  RDFWriteLayersClosest(&rdf, "rdfclo");
  RDFMolecWrite(&rdf, "rdfmol");
  RDFMolecWriteLayers(&rdf, "rdfmol");
  RDFMolecWriteLayersClosest(&rdf, "rdfmolclo");
  DoCWrite(&rdf, "DoC");
  DoCHistWrite(&rdf, "DoCHist");
  CoordNumWrite(&rdf, "coordnum");
  CoordNumHistWrite(&rdf, "coordnumHist");
  ACWriteAtomCounts(&ac, "atoms");
  ACWriteDensity(&ac, "density");
  ACWriteElectrodeCharge(&ac, "elecharge");
  ACWriteIons(&ac, "ions");
  ACWriteIonsInLayers(&ac, "layers");
  ACWriteIonsInLayersTime(&ac, "numionslayers");
  ACWriteElecChargeSlicesTime(&ac, "elecchargeslices");
  ACWriteElecNumIonsSlicesTime(&ac, "numionselecslices");
  ACWriteCollectiveVars(&ac, "ionCV");
  ACWriteElecAtomChargeHist(&ac, "elecatomchargehist");
  ACWriteElecAtomCharge(&ac, "elecatomcharge");
}
