#include <cstring>
#include "global.h"
#include "Reader.h"
#include "F1A_solver.h"

using namespace std;

int main(int argc, char*argv[])
{

  if (argc == 1)
  {
    Fatal_Error("Unspecified input file");
  }
  else
  {
    if (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help"))
    {
      cout << "FWH solver using Farassat 1A formulation" << endl;
      cout << "Developers: Trushant Patel, Weiqi Shen." << endl;
      cout << "Theoretical Fluid Dynamics and Turbulence Group" << endl;
      cout << "   Mechanical & Aerospace Engineering, UFL     " << endl;
      return 0;
    }
    else if (!strcmp(argv[1], "-v") || !strcmp(argv[1], "--version"))
    {
      cout << "Version: " << VERSION << endl;
      return 0;
    }
  }

  // create input file reader
  Reader input_reader(argv[1]);
  input_reader.read();

  // Finds out the pressure at the observer position by Farassat's 1A formulation.
  F1A_solver f1a;
  f1a.solve();

  return 0;
}
