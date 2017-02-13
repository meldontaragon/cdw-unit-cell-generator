/*
  Copyright (C) 2016-2017 David C. Miller

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

/*
  Input: 
  num - number of sites to create
  orig_lattice - original lattice of the high symmetry primitive unit cell
  supercell - super cell to use for new unit cell
  randomize - 1 to slightly randomize metal ion location and 0 to use high
  symmetry locations
  inversion - 1 to use inversion symmetry (1T) and 0 to use reflection symmetry
  (2H) 
  layers - number of layers (0 for bulk, 1 for monolayer, etc.)

  Output:
  atoms_m - Cartesian coordinates for metal ion sites (array of Location struct
  which holds an x, y, and z coordinate along with the atomic symbol) size num
 */
void make_m_site\
(Location atoms_m[], const unsigned num, double orig_lattice[3],\
 int supercell[2][2], const int randomize, const int inversion,\
 const unsigned layers)
{
  /* used for randomization of coordinates */
  double r1 = 0, r2 = 0;
  /* scale of randomization (default is 3%) */
  double min_scale = orig_lattice[0] * 0.01;
  double max_scale = orig_lattice[0] * 0.04;
  int sign = 1;

  unsigned ii;

   /* fractional location of atom within supercell */
  double (*frac_loc)[3];
  
  frac_loc = malloc(sizeof(double[3])*num);
  
  generate_frac_coord(frac_loc, num, supercell);

  srand((unsigned)(time(NULL)));
  for (ii = 0; ii < num; ii++)
    {
      if (randomize)
	{
	  /* 
	     generates a random value in [-b,-a] U [a,b] where a is a
	     minimum scaling and b is a maximum scaling
	  */

	  /* first we get a value between min_scale and max_scale */
	  r1 = ( (1.0*rand()/RAND_MAX)*(max_scale-min_scale) ) + min_scale;
	  r2 = ( (1.0*rand()/RAND_MAX)*(max_scale-min_scale) ) + min_scale;

	  /* gets a random value of either -1 or 1 to switch r1 and r2 */
	  sign = get_sign(((1.0*rand()/RAND_MAX)-0.5));
	  r1 = r1*sign;
	  sign = get_sign(((1.0*rand()/RAND_MAX)-0.5));
	  r2 = r2*sign;
	}
      
      /* 
	 location is n1*a1 + n2*a2 but a2=|a2|*0.5*x -sqrt(3)*0.5*|a2| 
	 r1 and r2 are the randomization values in the a and b directions 
      */

      if (!(inversion) && !(layers == 1))
	{
	  /* itr 0 */
	  atoms_m[ii].x = (frac_loc[ii][0] * orig_lattice[0])\
	    - (frac_loc[ii][1]*0.5*orig_lattice[0]) + r1;

	  atoms_m[ii].y = (sqrt(3.0)*0.5*frac_loc[ii][1] * orig_lattice[1]) +r2;

	  atoms_m[ii].z = 0.25 * orig_lattice[2];

	  /* itr 1 */
	  atoms_m[ii+num].x = (frac_loc[ii][0] * orig_lattice[0])\
	    - (frac_loc[ii][1]*0.5*orig_lattice[0]) + r1;

	  atoms_m[ii+num].y = (sqrt(3.0)*0.5*frac_loc[ii][1] * orig_lattice[1]) +r2; 

	  atoms_m[ii+num].z = 0.75 * orig_lattice[2];
	}
      else
	{
	  atoms_m[ii].x =  (frac_loc[ii][0] * orig_lattice[0])\
	    - (frac_loc[ii][1]*0.5*orig_lattice[0]) + r1;

	  atoms_m[ii].y = (sqrt(3.0)*0.5*frac_loc[ii][1] * orig_lattice[1]) +r2;

	  atoms_m[ii].z = 0.5 * orig_lattice[2];
	}
    }

  free (frac_loc);
}/* end of makeMSite(...) */

void make_x_site
(Location atomsX[], const unsigned num, double orig_lattice[3],\
 int supercell[2][2], const int inversion, const unsigned layers)
{
  /* correction vector for displacement of X atoms relative to M sites */
  double correction[2][2];

  unsigned ii;
  /* fractional location of atom within supercell */
  double (*frac_loc)[3];
  
  frac_loc = malloc(sizeof(double[3])*num);
  generate_frac_coord(frac_loc, num, supercell);
  
  for (ii = 0; ii < num; ii++)
    {
      /* generates the correction vectors for 1T and 2H phases */
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
	  /* itr = 0 */
	  atomsX[ii].x =	((frac_loc[ii][0] + correction[0][0]) * orig_lattice[0])\
	    - ((frac_loc[ii][1] + correction[0][1])*0.5*orig_lattice[0]);

	  atomsX[ii].y = (frac_loc[ii][1] + correction[0][1])\
	    *sqrt(3)*0.5*orig_lattice[1];

	  atomsX[ii].z = 0.125 * orig_lattice[2];
      
	  /* itr = 1 */
	  atomsX[ii+num].x = ((frac_loc[ii][0] + correction[1][0]) * orig_lattice[0]) \
	    - ((frac_loc[ii][1] + correction[1][1])*0.5*orig_lattice[0]);

	  atomsX[ii+num].y = (frac_loc[ii][1] + correction[1][1])\
	    *sqrt(3)*0.5*orig_lattice[1];

	  atomsX[ii+num].z = 0.375 * orig_lattice[2];
	  
	  /* itr = 2 */
	  atomsX[ii+2*num].x = ((frac_loc[ii][0] - correction[0][0]) * orig_lattice[0])\
	    - ((frac_loc[ii][1] - correction[0][1])*0.5*orig_lattice[0]);

	  atomsX[ii+2*num].y = (frac_loc[ii][1] - correction[0][1])\
	    *sqrt(3)*0.5*orig_lattice[1];

	  atomsX[ii+2*num].z = 0.625 * orig_lattice[2];
	    
	  /* itr = 3 */
	  atomsX[ii+3*num].x  = ((frac_loc[ii][0] - correction[1][0]) * orig_lattice[0])\
	    - ((frac_loc[ii][1] - correction[1][1])*0.5*orig_lattice[0]);

	  atomsX[ii+3*num].y = (frac_loc[ii][1] - correction[1][1])\
	    *sqrt(3)*0.5*orig_lattice[1];

	  atomsX[ii+3*num].z = 0.875 * orig_lattice[2];
	}
      else
	{
	  /* itr = 0 */
	  atomsX[ii].x = ((frac_loc[ii][0] + correction[0][0]) * orig_lattice[0])\
	    - ((frac_loc[ii][1] + correction[0][1])*0.5*orig_lattice[0]);
	  atomsX[ii].y = (frac_loc[ii][1] + correction[0][1]) * sqrt(3)*0.5*orig_lattice[1];
	  atomsX[ii].z = 0.75 * orig_lattice[2];
      
	  /* itr = 1 */
	  atomsX[ii+num].x = ((frac_loc[ii][0] + correction[1][0]) * orig_lattice[0]) \
	    - ((frac_loc[ii][1] + correction[1][1])*0.5*orig_lattice[0]);
	  atomsX[ii+num].y = (frac_loc[ii][1] + correction[1][1]) * sqrt(3)*0.5*orig_lattice[1];
	  atomsX[ii+num].z = 0.25 * orig_lattice[2];
	}
    
    }
  free (frac_loc);
}/* end of makeXSite(...) */
