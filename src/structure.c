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
 int strained, int strain_axis[3], int strain_min, int strain_max)
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
  double unstrained_lattice[3];
  
  char name[100], sym_type[20], inv_type[20], s_elM[5], s_elX[5], layer_type[40], filename[100], cdw_type[40];

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
  num = (unsigned) ceil(sqrt(supercell[0][0]*supercell[0][0] \
			     + supercell[0][1]*supercell[0][1] \
			     - supercell[0][0]*supercell[0][1])	\
			*sqrt(supercell[1][0]*supercell[1][0] \
			      + supercell[1][1]*supercell[1][1]	\
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

  unstrained_lattice[0] = orig_lattice[0];
  unstrained_lattice[1] = orig_lattice[1];
  unstrained_lattice[2] = orig_lattice[2];

  /* actual structure generation contained here */
  if (strained)
    {
      printf("Adding strain\n");
      strain_end = strain_max;
      strain_start = strain_min;
    }
  else
    {
      printf("Not strained\n");
      strain_start = 0;
      strain_end = 0;
    }
  
  for (i = strain_start; i <= strain_end; ++i)
    {
      printf("Current strain: %d percent\n",i);
      if (strain_axis[0])
	delta_x = (1.0) + i*0.01;
      if (strain_axis[1])
	delta_y = (1.0) + i*0.01;
      if (strain_axis[2])
	delta_z = (1.0) + i*0.01;

      /* printf("orig_lattice[0]: %f\n", unstrained_lattice[0]); */
      orig_lattice[0] = unstrained_lattice[0] * delta_x; 
      orig_lattice[1] = unstrained_lattice[1] * delta_y; 
      orig_lattice[2] = unstrained_lattice[2] * delta_z; 

      /* printf(" new_lattice[0]: %f\n\n\n", orig_lattice[0]); */
      /* make the M and X sites */
      makeMSite(atomsM, num, orig_lattice, supercell, randomize, inversion, layers);
      makeXSite(atomsX, num, orig_lattice, supercell, inversion, layers);

      /* get the new lattice parameters */
      lattice[0][0] = supercell[0][0]*unstrained_lattice[0]\
	- 0.5*supercell[0][1]*unstrained_lattice[0];
      lattice[0][1] = sqrt(3)*0.5*supercell[0][1]*unstrained_lattice[1];

      lattice[1][0] = supercell[1][0]*unstrained_lattice[0]\
	- 0.5*supercell[1][1]*unstrained_lattice[0];
      lattice[1][1] = sqrt(3)*0.5*supercell[1][1]*unstrained_lattice[1];

      lattice[0][0] *= delta_x;
      lattice[0][1] *= delta_x;

      lattice[1][0] *= delta_y;
      lattice[1][1] *= delta_y;

      if (layers == 1)
	{
	  lattice[2][2] = MONOLAYER_C_LAT;
	}
      else if (layers == 0)
	{
	  if (inversion)
	    {
	      lattice[2][2] = unstrained_lattice[2];
	      lattice[2][2] *= delta_z;
	    }
	  else
	    {
	      lattice[2][2] = unstrained_lattice[2];
	      lattice[2][2] *= delta_z;
	      num *= 2;
	    }
	}
      else
	{
	  exit(-1);
	}

      if (randomize) strcpy(sym_type,"Rand");
      else strcpy(sym_type, "HS");

      /*
	this program uses the naming convention that
	even monolayer 2H structures use 2H as opposed
	to 1H
      */
      if (inversion) strcpy(inv_type, "1T");
      else strcpy(inv_type, "2H");

      strcpy(s_elM, Symbol[elemM]);
      strcpy(s_elX, Symbol[elemX]);

      if (layers == 1) strcpy(layer_type,"Monolayer");
      else if (layers == 2) strcpy(layer_type,"Bilayer");
      else if (layers == 3) strcpy(layer_type,"Trilayer");
      else if (layers == 0) strcpy(layer_type,"Bulk");
      else strcpy(layer_type,"Other");

      /*
	the file name should match up with the standard I've started using:
	name = $INV_TYPE-$ELM$ELX2-$CDW_TYPE-Monolayer/Bulk/Bilayer-HS/Rand.vasp
	strain_name = Strain_PERC--$NAME
	
	IMPORTANT
	Here CDW_TYPE is *not* as easily defined from the input parameters
	i.e. (3,0),(0,3) is 3x3 while the Star of David CDW is more complex
	taking the form (4,1),(-1,3)

	My current (and poor) method is to hard code the cases I need in but
	this explicitly requires the code to be modified for any new cases.
	Any (a,0),(0,b) cases should work fine (since these are just axb)
	but special cases will result in a different file naming format.	
      */

      if ( (supercell[0][1] == 0) && (supercell[1][0] == 0) )
	sprintf(cdw_type,"%dx%d",supercell[0][0],supercell[1][1]);
      else if ( (supercell[0][0] == 4 ) && (supercell[0][1] == 1) && (supercell[1][0] == -1) && (supercell[1][1] == 3) )
	sprintf(cdw_type,"SoD");
      else
	sprintf(cdw_type,"(%d,%d)x(%d,%d)",supercell[0][0],supercell[0][1],supercell[1][0],supercell[1][1]);
	
      if (i != 0)
	sprintf(name,"strain_%+d--%s-%s%s2-%s %s %s",i,inv_type, s_elM, s_elX, layer_type, \
		sym_type, cdw_type);
      else
	sprintf(name,"%s-%s%s2-%s %s %s", inv_type, s_elM, s_elX, layer_type, \
		sym_type, cdw_type);
    
      if (i != 0)
	sprintf(filename,"Strain_%+d--%s-%s%s2-%s-%s-%s.vasp",i,inv_type,	\
		s_elM, s_elX, cdw_type, layer_type, sym_type);
      else
	sprintf(filename,"%s-%s%s2-%s-%s-%s.vasp",inv_type, s_elM, s_elX,\
		cdw_type, layer_type, sym_type);

      /* print to file filename */
      print_VASP_to_file(atomsM, atomsX, num, lattice, name, s_elM, s_elX, filename);
    }  
  /* clean up memory */
  free (atomsM);
  free (atomsX);
  
  return 1;
} /* int makeStructure(arrays...) */



