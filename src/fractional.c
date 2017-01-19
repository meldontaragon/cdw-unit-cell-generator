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

/* **************************************************************************
   fractional.c
   code written by David Miller
   mill2723 at msu dot edu

   Generates the fractional coordinates for a TMD trilayer. These coordinates
   are for the high symmetry locations of the metal ions (M in MX2 formula). For
   the chalcogen sites (X) the resulting coordinate needs to be modified. This
   is performed in the `make_m_site(...)` function.
   
*****************************************************************************/

#include "header.h"

/*
  Input:
  num - number of fractional coordinates to generate
  supercell - super cell to use for generation

  Output:
  frac_loc - 2d double array of size num by 3 with fractional
  coordinates
 */

void generate_frac_coord\
(double frac_loc[][3], const unsigned num, const int supercell[2][2])
{
  unsigned xmax = 0, ymax = 0, count = 0, ix, jy, kk;
  double min_angle;

  /* maximum x and y coordinates */
  xmax = (unsigned) ( abs(supercell[0][0]) + abs(supercell[1][0]) );
  ymax = (unsigned) ( abs(supercell[0][1]) + abs(supercell[1][1]) );

  /* minimum angle of vector for coordinates that lie in unit cell of output */
  min_angle = get_lattice_vector_angle(supercell[0][0], supercell[0][1]);

  for (ix = 0; ix < xmax; ++ix)
    {
      for (jy = 0; jy < ymax; ++jy)
	{
	  if (count >= num) { break; }	  
	  
	  if ( (ix == 0) && (jy == 0) )
	    {
	      /* first fraction coordinates are always at (0,0,0) */
	      for (kk = 0; kk < 3; ++kk)
		  frac_loc[0][kk] = 0;
	      /* increment the number of coordinates created */
	      count++;
	    }
	  else
	    {
	      /*
		checks if the angle of the new fractional location lies in the region
		where atoms are being added

		this is determined by a minimum angle for the vector (there is no
		distance requirement imposed, this is managed by xmax and ymax)
	      */
	      if( ( min_angle - get_lattice_vector_angle((int)ix, (int)jy) ) < EPS )
		    {
		      frac_loc[count][0] = ix; 
		      frac_loc[count][1] = jy; 
		      frac_loc[count][2] = 0;
  		    count++;
		    }
	    } /* end of else */  
	  } /* end of for on jy */
  } /* end of for on ix */
  if (count < num)
    fprintf(stderr,"Not enough coordinates made...\n");
} /* end of generate_frac_coord(...) */
