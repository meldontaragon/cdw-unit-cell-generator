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
   generateTMDSuperCell.c
   code written by David Miller
   mill2723 at msu dot edu

   Generates coordinates for arbitrarily sized supercells of TMDs (MX2)
   layered materials.
*********************************************************************/
#include "header.h"

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
    
    12 Strain T/F
    13 Strain axis a T/F
    14 Strain axis b T/F
    15 Strain axis c T/F
  */
  /* (1-2) */
  double orig_lattice[3];

  /* (3-6) */
  int supercell[2][2];

  /* (7) */
  /* generates bulk if 0, monolayer if 1, and nothing else is implemented */
  unsigned layers = 0;
  
  /* (8) */
  /* 1T has inversion symmetry, 1H has only reflection through xy plane */
  int inversion = 1;

  /* (9) */
  /* Choice of whether to randomize the coordinates for CDW searches */
  int randomize = 0;

  /* (10-11) */
  AtomicSymbol elemM = Ta, elemX = S;

  /* (12) */
  int strained = 0;

  /* (13-15) */
  int strain_axis[3] = {0, 0, 0};
  if (argc <= 11)
    {
      printf("Not the correct number of parameters\n");
      printHelp();
      return -1;
    }  

  if (argc > 11)
    {
      orig_lattice[0] = atof(argv[1]);
      orig_lattice[1] = orig_lattice[0];
      orig_lattice[2] = atof(argv[2]);

      supercell[0][0] = atoi(argv[3]);
      supercell[0][1] = atoi(argv[4]);
      supercell[1][0] = atoi(argv[5]);
      supercell[1][1] = atoi(argv[6]);

      layers = (unsigned) (atoi(argv[7]));

      inversion = atob(argv[8][0]);
      randomize = atob(argv[9][0]);

      elemM = (AtomicSymbol) (atoi(argv[10])-1);
      elemX = (AtomicSymbol) (atoi(argv[11])-1);
	}
  if (argc > 15)
    {
      strained = atob(argv[12][0]);
      strain_axis[0] = atob(argv[13][0]);
      strain_axis[1] = atob(argv[14][0]);
      strain_axis[2] = atob(argv[15][0]);
    }
   
  makeStructure(orig_lattice, supercell, inversion, randomize, layers,\
		elemM, elemX, strained, strain_axis);
  return 1;
}

