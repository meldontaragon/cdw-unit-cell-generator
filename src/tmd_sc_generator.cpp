/*
  This file is part of TMD CDW Unit Cell Generator

  TMD CDW Unit Cell Generatort is free software: you can redistribute it
  and/or modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation, either version
  3 of the License, or (at your option) any later version.

  TMD CDW Unit Cell Generator is distributed in the hope that it will be
  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with TMD CDW Unit Cell Generator.  If not, see
  <http://www.gnu.org/licenses/>.
*/
/* ******************************************************************
   tmd_sc_generator.cpp
   code written by David Miller
   mill2723 at msu dot edu


   Generates coordinates for arbitrarily sized supercells of MX2
   layered materials.
*********************************************************************/
#ifndef MONOLAYER_C_LAT
#define MONOLAYER_C_LAT 23.6000
#endif

#include "header.h"

using namespace std;

int main(int argc, char* argv[])
{
  /*
    read order:
    1  Orig lattice a
    2  Orig lattice c
    3  Supercell a1
    4  Supercell a2
    5  Supercell b1
    6  Supercell b2
    7  Monolayer T/F
    8  Inversion T/F
    9  Randomize T/F
    10 Elements M (Ta,Ti) ! use atomic number
    11 Elements X (S,Se,Te) ! use atomic number
  */
  vector<double> orig_lattice(3,0.0);
  vector< vector<int> > supercell(2,vector<int>(2,0));

  /* default values */
  
  /* 1T has inversion symmetry, 1H has only reflection through xy plane */
  bool inversion = true;
  /* Choice of whether to randomize the coordinates for CDW searches */
  bool randomize = false; 
  bool monolayer = true;
  Elements::AtomicSymbol elemM = Elements::Ta, elemX = Elements::S;

  if (argc != 12)
    {
      cerr << "Not the correct number of parameters" << endl;
      printHelp();
      return -1;
    }
  else
    {
      orig_lattice[0] = atof(argv[1]);
      orig_lattice[1] = orig_lattice[0];
      orig_lattice[2] = atof(argv[2]);

      supercell[0][0] = atoi(argv[3]);
      supercell[0][1] = atoi(argv[4]);
      supercell[1][0] = atoi(argv[5]);
      supercell[1][1] = atoi(argv[6]);

      switch(argv[7][0])
	{
	case 'T':
	case 't':
	  monolayer = true;
	  break;
	case 'F':
	case 'f':
	  monolayer =false;
	  break;
	default:
	  cerr << "Incorrect logical statement." << endl;
	  return -1;
	}

      switch(argv[8][0])
	{
	case 'T':
	case 't':
	  inversion = true;
	  break;
	case 'F':
	case 'f':
	  inversion =false;
	  break;
	default:
	  cerr << "Incorrect logical statement." << endl;
	  return -1;
	}

      switch(argv[9][0])
	{
	case 'T':
	case 't':
	  randomize = true;
	  break;
	case 'F':
	case 'f':
	  randomize =false;
	  break;
	default:
	  cerr << "Incorrect logical statement." << endl;
	  return -1;
	}

      elemM = Elements::AtomicSymbol(atoi(argv[10])-1);
      elemX = Elements::AtomicSymbol(atoi(argv[11])-1);
    }
  
  makeStructure(orig_lattice, supercell, inversion, randomize, monolayer, elemM, elemX);
}

int makeStructure(vector<double> & orig_lattice, vector< vector<int> > & supercell,
		  bool inversion, bool randomize, bool monolayer,
		  Elements::AtomicSymbol elemM, Elements::AtomicSymbol elemX)
{
  /* ******************************************
     NEEDED INFORMATION AND FORMAT
   
     (A) Original unit cell
     - Hexagonal unit cell requires
     - a and c lattice parameters
     (B) 1T vs 1H vs 2H etc. 
     - currently only monolayers
     (C) Choice of atoms for M and X site
     - defaults for TaS2
     (D) Choice of supercell 
     - given in {{xa,xb},{ya,yb}} format
     (E) High symmetry vs random pertubations
     - necessary for CDW searchs that
     - require length contractions
  ****************************************/

  /* number of M atoms, number of X atoms = 2*num */
  unsigned num = 0; 
  
  vector<Location> atomsM = vector<Location>(num);
  vector<Location> atomsX = vector<Location>(2*num);

  /* currently always uses a hexagonal primitive cell */
  vector<double> angle(3);
   
  angle[0] = 90.;
  angle[1] = angle[0];
  angle[2] = 120.;

  /*
    num = sqrt(a' * a') * sqrt(b' * b') / sqrt(a*a)
    where the a' and b' vectors are the supervectors
    defined by (xa, xb), (ya, yb) where a' = xa*a + xb*b
    and b' = ya*a + yb*b and sqrt(a*a)=sqrt(b*b)
  */

  /*
    since a*b = (1/2)*|a|*2 for hexagonal unit cells the 
    (x+y)^2 = x^2+y^2+2xy ---> x^2+y^2+|x||y| 
  */
  num = ceil(sqrt(supercell[0][0]*supercell[0][0]
		  + supercell[0][1]*supercell[0][1]
		  - supercell[0][0]*supercell[0][1])
	     *sqrt(supercell[1][0]*supercell[1][0]
		   + supercell[1][1]*supercell[1][1]
		   - supercell[1][0]*supercell[1][1]));

  /* debug statement
     cout << "a'^2: " << (supercell[0][0]*supercell[0][0] + \
     supercell[0][1]*supercell[0][1] - supercell[0][0]*supercell[0][1]) \
     << "\t b'^2: " << (supercell[1][0]*supercell[1][0] + supercell[1][1]*supercell[1][1]\
     - supercell[1][0]*supercell[1][1]) << endl;
  */
    
  /* debug statement
     cout << "DEBUG: Making M site locations..." << endl;
  */
  makeMSite(atomsM, num, orig_lattice, supercell, randomize, inversion, monolayer);
  /* debug statement
     cout << "DEBUG: Making X site locations..." << endl;
  */
  makeXSite(atomsX, num, orig_lattice, supercell, inversion, monolayer);

  vector< vector<double> > lattice = vector< vector<double> >(3,vector<double>(3,0.0));
  /*
    need to generate the new lattice for output
    angles may also change now (although I don't believe they will)
  */
  lattice[0][0] = supercell[0][0]*orig_lattice[0]
    - 0.5*supercell[0][1]*orig_lattice[0];
  lattice[0][1] = sqrt(3)*0.5*supercell[0][1]*orig_lattice[1];

  lattice[1][0] = supercell[1][0]*orig_lattice[0]
    - 0.5*supercell[1][1]*orig_lattice[0];
  lattice[1][1] = sqrt(3)*0.5*supercell[1][1]*orig_lattice[1];

  if (monolayer)
    {
      lattice[2][2] = MONOLAYER_C_LAT;
    }
  else
    {
      if (inversion)
	{
	  lattice[2][2] = orig_lattice[2];
	}
      else
	{
	  lattice[2][2] = orig_lattice[2];
	  num *= 2;
	}
    }
    
  string name;
  string sym_type, inv_type, s_elM, s_elX, layer_type;
  if (randomize) sym_type = "Randomized";
  else sym_type = "High Symmetry";

  if (inversion) inv_type = "1T";
  else inv_type = "2H";

  s_elM = Elements::Symbol[elemM];
  s_elX = Elements::Symbol[elemX];

  if (monolayer) layer_type = "monolayer";
  else layer_type = "bulk";

  char buffer[200];
  sprintf(buffer,"%s-%s%s2-%s %s (%d,%d)x(%d,%d)", inv_type.c_str(),
	  s_elM.c_str(),s_elX.c_str(),layer_type.c_str(),sym_type.c_str(),
	  supercell[0][0],supercell[0][1],supercell[1][0],supercell[1][1]);
  name = string(buffer);
  printVASP(atomsM, atomsX, num, lattice, name, s_elM, s_elX);

  return 1;
} //int main(...)

