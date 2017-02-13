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

/* make sites */
void make_m_site(Location atoms_m[], const unsigned num, double orig_lattice[3],\
		 int supercell[2][2], const int randomize, const int inversion,\
		 const unsigned layers);
void make_x_site(Location atomsX[], const unsigned num, double orig_lattice[3],\
		 int supercell[2][2], const int inversion, const unsigned layers);

/* output */
void print_vasp_to_file(Location loc_m[], Location loc_x[], unsigned n,\
			double lattice[3][3], char * name, char * elem_m,\
			char * elem_x, char * file_name);

void print_xyz_to_file(Location loc_m[], Location loc_x[], unsigned n,\
		       char * filename) __attribute__((deprecated));

/* these functions should not be used anymore */
void print_vasp(Location loc_m[], Location loc_x[], unsigned n,\
	       double lattice[3][3], char * name, char * elem_m, char * elem_x)\
  __attribute__((deprecated)); 
void print_xyz(Location loc_m[], Location loc_x[], unsigned n)\
  __attribute__((deprecated));

void print_help();
void print_version();
void print_test_start();

/* fractional coordinate generation */
void generate_frac_coord(double frac_loc[][3], const unsigned num,\
			 int supercell[2][2]);

/* structure generation functions */
int make_structure(double orig_lattice[3], int supercell[2][2],\
		   const int inversion, const int randomize, const unsigned layers, \
		   const AtomicSymbol elem_m, const AtomicSymbol elem_x,\
		   const int strained, int strain_axis[3], const int strain_min, \
		   const int strain_max);

/* additional functions used throughout the code */
double get_lattice_vector_angle(const int a, const int b);

double dtor(double deg) __attribute__ ((const));
int atob(char a) __attribute__ ((pure));

int get_sign(const double value) __attribute__((const));

/* test suite related functions */
void test_print_vasp_monolayer();
void test_print_vasp_bulk();
void test_print_xyz_monolayer();
void test_print_xyz_bulk();
 
#endif
