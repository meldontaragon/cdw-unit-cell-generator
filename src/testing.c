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

/* ******************************************************************
   testing.c
   code written by David Miller
   mill2723 at msu dot edu

   Code to interface with test suite. This performs simple outputs
   which can be verified against expected results to check if the program is
   working as expected and desired
*********************************************************************/
#include "header.h"
void __test_print_vasp(const unsigned layers);
void __test_print_vasp_strain(const unsigned type);

void __test_print_vasp(const unsigned layers)
{
  double orig_lattice[3] = {5, 5, 10};
  int supercell[2][2] = { {1,0}, {0,1} };
  const AtomicSymbol elem_m = Ta, elem_x = S;

  const int strain = 0;
  int strain_logic[3] = {0,0,0};
  int abs_strain = 1;
  double strain_val = 0.01;

  /* test output for strain with -1% and 1% along with default unstrained file*/
  make_structure(orig_lattice, supercell, 1, 0, layers, elem_m, elem_x,\
		 strain, strain_logic, abs_strain, strain_val, 0, 0);

  make_structure(orig_lattice, supercell, 0, 0, layers, elem_m, elem_x,\
		 strain, strain_logic, abs_strain, strain_val, 0, 0);
}

void __test_print_vasp_strain(const unsigned type)
{
  /* type is 0: biaxial, 1: uniaxial */

  double orig_lattice[3] = {5, 5, 10};
  int supercell[2][2] = { {1,0}, {0,1} };
  const AtomicSymbol elem_m = Ta, elem_x = S;
  unsigned layers = 1;
  
  const int strain = 1;
  int strain_logic[3] = {0,0,0};
  int abs_strain = 1;
  double strain_val = 0.01;

  if (type == 0)
    {
      strain_logic[0] = 1;
      strain_logic[1] = 1;
    }
  else if (type == 1)
    {
      strain_logic[1] = 1;
    }

    /* test output for strain with -1% and 1% along with default unstrained file*/
  make_structure(orig_lattice, supercell, 1, 0, layers, elem_m, elem_x,\
		 strain, strain_logic, abs_strain, strain_val, -2, 2);

  make_structure(orig_lattice, supercell, 0, 0, layers, elem_m, elem_x,\
		 strain, strain_logic, abs_strain, strain_val, -2, 2);

}

void test_print_vasp_monolayer()
{
  __test_print_vasp((unsigned) 1);
}

void test_print_vasp_bulk()
{
  __test_print_vasp((unsigned) 0);
}

void test_print_vasp_biaxial_strain()
{
  __test_print_vasp_strain((unsigned) 0);
}

void test_print_vasp_uniaxial_strain()
{
  __test_print_vasp_strain((unsigned) 1);
}

void test_print_xyz_monolayer()
{
  printf("Not implemented yet...\n");
}

void test_print_xyz_bulk()
{
  printf("Not implemented yet...\n");
}
