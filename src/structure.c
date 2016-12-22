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

/* ************************************************************
   structure.c
   code written by David Miller
   mill2723 at msu dot edu

   Given a set of inputs this generates and outputs a
   set of structures (with strain added in) to files
   defined by the input parameters)
***************************************************************/
#include "header.h"

#ifndef MONOLAYER_C_LAT
#define MONOLAYER_C_LAT 23.6000
#endif

int makeStructure
(double orig_lattice[3], int supercell[2][2], int inversion, int randomize,\
 unsigned layers, AtomicSymbol elemM, AtomicSymbol elemX,\
 int strained, int strain_axis[3], int strain_min_max[2])
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
     (G) Amount of strain to put on the system
     - with a minimum and maximum percentage
  ****************************************/

  /* number of M atoms, number of X atoms = 2*num */
  unsigned num = 0;
  Location *atomsM, *atomsX;

  double lattice[3][3];
  double angle[3];
  int strain_start, strain_end, i;
  double delta_x = 1.0, delta_y = 1.0, delta_z = 1.0;
  
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

  /* actual structure generation contained here */
  if (strained)
    {
      strain_end = strain_min_max[1];
      strain_start = strain_min_max[0];
    }
  else
    {
      strain_start = 0;
      strain_end = 0;
    }
  
  for (i = strain_start; i <= strain_end; ++i)
    {
      if (strain_axis[0])
	delta_x = (1.0) + i*0.01;
      if (strain_axis[1])
	delta_y = (1.0) + i*0.01;
      if (strain_axis[2])
	delta_z = (1.0) + i*0.01;

      orig_lattice[0] = orig_lattice[0] * delta_x;
      orig_lattice[1] = orig_lattice[1] * delta_y;
      orig_lattice[2] = orig_lattice[2] * delta_z;
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

      if (randomize) strcpy(sym_type,"rand");
      else strcpy(sym_type, "hsym");

      if (inversion) strcpy(inv_type, "1T");
      else strcpy(inv_type, "2H");

      strcpy(s_elM, Symbol[elemM]);
      strcpy(s_elX, Symbol[elemX]);

      if (layers == 1) strcpy(layer_type,"monolayer");
      else if (layers == 2) strcpy(layer_type,"bilayer");
      else if (layers == 0) strcpy(layer_type, "bulk");
      else strcpy(layer_type,"unknown");

      if (i != 0)
	sprintf(buffer,"strain_%d--%s-%s%s2-%s %s (%d,%d)x(%d,%d)",i,inv_type, s_elM, s_elX, layer_type, \
		sym_type, supercell[0][0],supercell[0][1],supercell[1][0],supercell[1][1]);
      else
	sprintf(buffer,"%s-%s%s2-%s %s (%d,%d)x(%d,%d)", inv_type, s_elM, s_elX, layer_type, \
		sym_type, supercell[0][0],supercell[0][1],supercell[1][0],supercell[1][1]);
      strcpy(name, buffer);
    
      if (i != 0)
	sprintf(buffer,"strain_%d--%s-%s%s2-%s-%s_%d-%d_%d-%d.vasp",i,inv_type,	\
		s_elM, s_elX, layer_type, sym_type, supercell[0][0], supercell[0][1], supercell[1][0], supercell[1][1]);
      else
	sprintf(buffer,"%s-%s%s2-%s-%s_%d-%d_%d-%d.vasp",inv_type, s_elM, s_elX,\
		layer_type, sym_type, supercell[0][0], supercell[0][1], supercell[1][0], supercell[1][1]);

      print_VASP_to_file(atomsM, atomsX, num, lattice, name, s_elM, s_elX, buffer);
    }  
  /* clean up memory */
  free (atomsM);
  free (atomsX);
  
  return 1;
} /* int makeStructure(arrays...) */



