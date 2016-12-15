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
/* ********************************************************************
   atomic_sites.c
   code written by David Miller
   mill2723 at msu dot edu

   Generates atomic sites for a given lattice using fractional coordinates
   and sets them on a super lattice as defined by a set of vectors. This
   accounts for different TMD structures (1T and 2H) as well as slightly
   randomizing the positions of the metal if so desired.

***********************************************************************/

#include "header.h"

void makeMSite\
(Location atomsM[], unsigned num, double orig_lattice[3],\
 int supercell[2][2], int randomize, int inversion, unsigned layers)
{
  double r1 = 0, r2 = 0;
  double scale = orig_lattice[0] * 0.03;

  unsigned n;

   /* fractional location of atom within supercell */
  double (*frac_loc)[3];
  
  frac_loc = malloc(sizeof(double[3])*num);
  
  generateFracCoord(frac_loc, num, supercell);

  srand((unsigned)(time(NULL)));
  for (n = 0; n < num; n++)
    {
      if (randomize)
	{
	  r1 = (((1.0*rand()/RAND_MAX)*2.0)-1.0)*scale;
	  r2 = (((1.0*rand()/RAND_MAX)*2.0)-1.0)*scale;
	}
      
      /* location is n1*a1 + n2*a2 but a2=|a2|*0.5*x -sqrt(3)*0.5*|a2| */

      if (!(inversion) && !(layers == 1))
	{
	  atomsM[n].x = (frac_loc[n][0] * orig_lattice[0])\
	    - (frac_loc[n][1]*0.5*orig_lattice[0]) + r1;
	  atomsM[n].y = (sqrt(3)*0.5*frac_loc[n][1] * orig_lattice[1]) +r2;
	  atomsM[n].z = 0.25 * orig_lattice[2];

	  atomsM[n+num].x = (frac_loc[n][0] * orig_lattice[0])\
	    - (frac_loc[n][1]*0.5*orig_lattice[0]) + r1;
	  atomsM[n+num].y = (sqrt(3)*0.5*frac_loc[n][1] * orig_lattice[1]) +r2; 
	  atomsM[n+num].z =	      0.75 * orig_lattice[2];
	}
      else
	{
	  atomsM[n].x =  (frac_loc[n][0] * orig_lattice[0])		\
	    - (frac_loc[n][1]*0.5*orig_lattice[0]) + r1;
	  atomsM[n].y = (sqrt(3)*0.5*frac_loc[n][1] * orig_lattice[1]) +r2;
	  atomsM[n].z = 0.5 * orig_lattice[2];
	}
    }
}/* end of makeMSite(...) */

void makeXSite
(Location atomsX[], unsigned num, double orig_lattice[3],\
 int supercell[2][2], int inversion, unsigned layers)
{
  /* correction vector for displacement of X atoms relative to M sites */
  double correction[2][2];

  unsigned n;
  /* fractional location of atom within supercell */
  double (*frac_loc)[3];
  
  frac_loc = malloc(sizeof(double[3])*num);
  /*
  frac_loc = (double**) malloc(sizeof(double*)*num);
  for (i = 0; i < num; ++i)
    frac_loc[i] = (double*) malloc(sizeof(double)*3);
  */
  generateFracCoord(frac_loc, num, supercell);
  
  for (n = 0; n < num; n++)
    {
      if (inversion)
	{
	  correction[0][0] = 1.0/3.0;
	  correction[0][1] = 2.0/3.0;
	  correction[1][0] = -1.0/3.0;
	  correction[1][1] = -2.0/3.0;
	}
      else
	{
	  correction[0][0] = 1.0/3.0;
	  correction[0][1] = 2.0/3.0;
	  correction[1][0] = 1.0/3.0;
	  correction[1][1] = 2.0/3.0;
	}
      
      if (!(inversion) && !(layers == 1))
	{
	  atomsX[n].x =	((frac_loc[n][0] + correction[0][0]) * orig_lattice[0])\
	    - ((frac_loc[n][1] + correction[0][1])*0.5*orig_lattice[0]);
	  atomsX[n].y = (frac_loc[n][1] + correction[0][1]) * sqrt(3)*0.5*orig_lattice[1];
	  atomsX[n].z = 0.125 * orig_lattice[2];
      
	  atomsX[n+num].x = ((frac_loc[n][0] + correction[1][0]) * orig_lattice[0]) \
	    - ((frac_loc[n][1] + correction[1][1])*0.5*orig_lattice[0]);
	  atomsX[n+num].x = (frac_loc[n][1] + correction[1][1]) * sqrt(3)*0.5*orig_lattice[1];
	  atomsX[n+num].x = 0.375 * orig_lattice[2];
	  
	  atomsX[n+2*num].x = ((frac_loc[n][0] - correction[0][0]) * orig_lattice[0])\
	    - ((frac_loc[n][1] - correction[0][1])*0.5*orig_lattice[0]);
	  atomsX[n+2*num].y = (frac_loc[n][1] - correction[0][1]) * sqrt(3)*0.5*orig_lattice[1];
	  atomsX[n+2*num].z = 0.625 * orig_lattice[2];
	    
	  atomsX[n+3*num].x  = ((frac_loc[n][0] - correction[1][0]) * orig_lattice[0]) \
	    - ((frac_loc[n][1] - correction[1][1])*0.5*orig_lattice[0]);
	  atomsX[n+3*num].y = (frac_loc[n][1] - correction[1][1]) * sqrt(3)*0.5*orig_lattice[1];
	  atomsX[n+3*num].z = 0.875 * orig_lattice[2];
	}
      else
	{
	  atomsX[n].x = ((frac_loc[n][0] + correction[0][0]) * orig_lattice[0])\
	    - ((frac_loc[n][1] + correction[0][1])*0.5*orig_lattice[0]);
	  atomsX[n].y = (frac_loc[n][1] + correction[0][1]) * sqrt(3)*0.5*orig_lattice[1];
	  atomsX[n].z = 0.75 * orig_lattice[2];
      
	  atomsX[n+num].x = ((frac_loc[n][0] + correction[1][0]) * orig_lattice[0]) \
	    - ((frac_loc[n][1] + correction[1][1])*0.5*orig_lattice[0]);
	  atomsX[n+num].y = (frac_loc[n][1] + correction[1][1]) * sqrt(3)*0.5*orig_lattice[1];
	  atomsX[n+num].z = 0.25 * orig_lattice[2];
	}
    
    }
  free (frac_loc);
}/* end of makeXSite(...) */
