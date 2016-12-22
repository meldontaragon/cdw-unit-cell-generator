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

/* ************************************************************************
   header_c.h
   code written by David Miller
   mill2723 at msu dot edu

   Contains all includes and definitions for additional files used
   throughout this code.

***************************************************************************/

#ifndef HEADER_C_H
#define HEADER_C_H
  
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "elements.h"

struct Location_S
{
  double x, y, z;
  AtomicSymbol elem;
};

typedef struct Location_S Location;

/* defines pi for future use */
static const double PI = 3.1415926535;
/* defines a minimum value for double comparision */
static const double EPS = 10E-5;

/* need to add __attribute__ ((deprecated)) to applicable functions */

/* make sites */
void makeMSite(Location atomsM[], unsigned num, double orig_lattice[3],\
	       int supercell[2][2], int randomize, int inversion, unsigned layers);
void makeXSite(Location atomsX[], unsigned num, double orig_lattice[3],\
	       int supercell[2][2], int inversion, unsigned layers);

/* output */
void printVASP(Location locTa[], Location locS[], unsigned n,\
	       double lattice[3][3], char * name, char * elemM, char * elemX);
void print_VASP_to_file(Location locTa[], Location locS[], unsigned n,\
			double lattice[3][3], char * name, char * elemM,\
			char * elemX, char * file_name);

void printXYZ(Location LocTa[], Location locS[], unsigned n);
void printHelp();

/* fractional coordinate generation */
void generateFracCoord(double frac_loc[][3], unsigned num,\
		       int supercell[2][2]);

/* structure generation functions */
int makeStructure(double orig_lattice[3], int supercell[2][2],\
		  int inversion, int randomize, unsigned layers,\
		  AtomicSymbol elemM, AtomicSymbol elemX, int strained,\
		  int strain_axis[3]);

double getLatticeVectorAngle(const int a, const int b);

double dtor(double deg) __attribute__ ((const));
int atob(char a) __attribute__ ((pure));
 
#endif
