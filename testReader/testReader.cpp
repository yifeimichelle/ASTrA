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
    if(argc != 2)
    {
        cout << "this program takes one argument that is the input file ";
        cout << endl;
        return 1;
    }

    cout << "should print \"comment found!\":" << endl;
    char firstChar[] = "#";
    char ignoreChar[] = "#";
    if( strcmp(firstChar,ignoreChar) != 0 )
      {
	cout << "not a comment." << endl;
      }
    else
      {
	cout << "comment found!" << endl;
      }
    cout << endl;
    
    // read inputs
    cout << "reading inputs..." << endl;
    string inputFile(argv[1]);
    System system = System(inputFile); // declare system
 
    cout << system.getTrajFile() << endl;
    //cout << m_inputs[2][2] << endl;
    cout << system.getBoxDim(2) << endl;
    cout << system.getNumMolecsOfType(4) << endl;

}
