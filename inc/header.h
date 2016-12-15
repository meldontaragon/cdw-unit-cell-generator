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
   header.h
   code written by David Miller
   mill2723 at msu dot edu

   Contains all includes and definitions for additional files used
   throughout this code.

***************************************************************************/
  
#ifndef HEADER_H
#define HEADER_H

#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>

#include "elements.h"
#include "Location.h"

using namespace std;

/* defines pi for future use */
const double PI = 4*atan(1.0);
/* defines a minimum value for double comparision */
const double EPS = 10E-5;

/* need to add __attribute__ ((deprecated)) to applicable functions */

/* make sites */
void makeMSite(vector<Location>& atomsM, unsigned num, vector<double> & orig_lattice, \
	       vector< vector<int> > & supercell, bool randomize, bool inversion, bool monolayer);
void makeXSite(vector<Location>& atomsX, unsigned num, vector<double> & orig_lattice,\
	       vector< vector<int> > & supercell, bool inversion, bool monolayer);

/* output */
void printVASP(vector<Location>& locTa, vector<Location> & locS, unsigned n,\
	       vector< vector<double> > & lattice, string name);
void printVASP(vector<Location> & locTa, vector<Location> & locS, unsigned n,\
	       vector< vector<double> > & lattice, string name, string elemM, string elemX);

void printXYZ(vector<Location>& LocTa, vector<Location>& locS, unsigned n);

/* fractional coordinate generation */
void generateFracCoord(vector< vector<double> > &frac_loc, unsigned num,\
		       vector< vector<int> > &supercell);

/* structure generation functions */

int makeStructure(vector<double> & orig_lattice, vector< vector<int> > & supercell,\
		  bool inversion, bool randomize, bool monolayer, Elements::AtomicSymbol elemM,\
		  Elements::AtomicSymbol elemX);

double getLatticeVectorAngle(int a, int b);
void printHelp();
double dtor(double deg);
bool atob(char a);
 
#endif
