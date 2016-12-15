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
   tmd_sc_generator.c
   code written by David Miller
   mill2723 at msu dot edu


   Generates coordinates for arbitrarily sized supercells of MX2
   layered materials.
*********************************************************************/
#ifndef MONOLAYER_C_LAT
#define MONOLAYER_C_LAT 23.6000
#endif

#include "header_c.h"

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
  unsigned layers;
  
  /* (8) */
  /* 1T has inversion symmetry, 1H has only reflection through xy plane */
  int inversion;

  /* (9) */
  /* Choice of whether to randomize the coordinates for CDW searches */
  int randomize;

  /* (10-11) */
  AtomicSymbol elemM = Ta, elemX = S;

  /* (12) */
  int strained;

  /* (13-15) */
  int strain_axis[3];
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

      elemM = atoi(argv[10])-1;
      elemX = atoi(argv[11])-1;
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

int makeStructure
(double orig_lattice[3], int supercell[2][2], int inversion, int randomize,\
 unsigned layers, AtomicSymbol elemM, AtomicSymbol elemX,\
 int strained, int strain_axis[3])
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
     (F) Strain and whether it is applied to
     - a given axis
  ****************************************/

  /* number of M atoms, number of X atoms = 2*num */
  unsigned num = 0;
  Location *atomsM, *atomsX;

  double lattice[3][3];
  double angle[3];

  char name[100], sym_type[20], inv_type[20], s_elM[5], s_elX[5], layer_type[40], buffer[200];

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
  num = (unsigned) ceil(sqrt(supercell[0][0]*supercell[0][0]	\
			     + supercell[0][1]*supercell[0][1]		\
			     - supercell[0][0]*supercell[0][1])		\
			*sqrt(supercell[1][0]*supercell[1][0]		\
			      + supercell[1][1]*supercell[1][1]		\
			      - supercell[1][0]*supercell[1][1]));

/* if 2H is used for a bilayer or bulk system then num is effectively doubled */
if (!(inversion) && !(layers == 1))
  {
    atomsM = malloc(sizeof(Location)*num*2);
    atomsX = malloc(sizeof(Location)*num*4);
  }
 else
   {
     atomsM = malloc(sizeof(Location)*num);
     atomsX = malloc(sizeof(Location)*num*2);
   }

/* make the M and X sites */
makeMSite(atomsM, num, orig_lattice, supercell, randomize, inversion, layers);

makeXSite(atomsX, num, orig_lattice, supercell, inversion, layers);

/* get the new lattice parameters */
lattice[0][0] = supercell[0][0]*orig_lattice[0]\
  - 0.5*supercell[0][1]*orig_lattice[0];
lattice[0][1] = sqrt(3)*0.5*supercell[0][1]*orig_lattice[1];

lattice[1][0] = supercell[1][0]*orig_lattice[0]\
  - 0.5*supercell[1][1]*orig_lattice[0];
lattice[1][1] = sqrt(3)*0.5*supercell[1][1]*orig_lattice[1];

if (layers == 1)
  {
    lattice[2][2] = MONOLAYER_C_LAT;
  }
 else if (layers == 0)
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
 else
   {
     exit(-1);
   }
    
 if (randomize) strcpy(sym_type,"Randomized");
 else strcpy(sym_type, "High Symmetry");

 if (inversion) strcpy(inv_type, "1T");
 else strcpy(inv_type, "2H");

 strcpy(s_elM, Symbol[elemM]);
 strcpy(s_elX, Symbol[elemX]);

 if (layers == 1) strcpy(layer_type,"monolayer");
 else if (layers == 2) strcpy(layer_type,"bilayer");
 else if (layers == 0) strcpy(layer_type, "bulk");
 else strcpy(layer_type,"unknown");

 sprintf(buffer,"%s-%s%s2-%s %s (%d,%d)x(%d,%d)", inv_type, s_elM, s_elX, layer_type, \
	 sym_type, supercell[0][0],supercell[0][1],supercell[1][0],supercell[1][1]);
 strcpy(name, buffer);
 printVASP(atomsM, atomsX, num, lattice, name, s_elM, s_elX);

 /* clean up memory */
 free (atomsM);
 free (atomsX);
  
 return 1;
} /* int makeStructure(arrays...) */


