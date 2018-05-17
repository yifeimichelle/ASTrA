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
    //string chargeFile(argv[2]); //!!
    System system(inputFile); // initialize system
    //System system(inputFile,chargeFile); //!!
    RDF rdf(system); // initialize rdf
    AtomCounter ac(system); // initialize atomcounter
    Frame frame(system); // initialize trajectory frame reader

  
    cout << "initial charges" << endl;
  
    cout << "Reading trajectory ..." << endl;

    for (int i=0; i<10; i++)
      {
	frame.readStep(5); //!!
	frame.readCharges(); //!!
	cout << frame.getAtom(4610).getCharge() << endl;
	cout << frame.sumCharges(system, system.getAnodeID()) << endl;
	cout << frame.sumCharges(system, system.getCathodeID()) << endl;
	
      }
}
