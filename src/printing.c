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
   printing.c
   code written by David Miller
   mill2723 at msu dot edu

   Contains all functions for output including the VASP and xyz format
   as well as the help printed to terminal.
***********************************************************************/

/* VERSION DO NOT EDIT */
#ifndef VERSION
#define VERSION "0.0.0-error"
#endif
/* VERSION DO NOT EDIT */

#include "header.h"
void print_xyz(Location loc_m[], Location loc_x[], unsigned n)
{
  unsigned j;
 fprintf(stdout,"%d\n",3*n);
 fprintf(stdout,"c-cdw-TaS2 super cell SQRT(13)xSQRT(13)x1 (in angstroms)\n");

  for (j = 0; j < n; j++)
    {
     fprintf(stdout,"Ta        %.4f        %.4f        %.4f\n",\
	     loc_m[j].x, loc_m[j].y,loc_m[j].z);
     fprintf(stdout,"S        %.4f        %.4f        %.4f\n",\
	     loc_x[j].x, loc_x[j].y,loc_x[j].z);
     fprintf(stdout,"S        %.4f        %.4f        %.4f\n",\
	     loc_x[j+n].x, loc_x[j+n].y,loc_x[j+n].z);
    }
}

void print_vasp_to_file(Location loc_m[], Location loc_x[], unsigned n,\
		   double lattice[3][3], char * name, char * elemM,\
		   char * elemX, char * file_name) 
{
  FILE * fp;
  unsigned i;

  fp = fopen(file_name, "w");

  fprintf(fp,"%s\n",name);
  fprintf(fp,"%.1f\n",1.0);

  for (i = 0; i < 3; i++)
    fprintf(fp,"   %.9f                %.9f                %.9f\n",	\
	   lattice[i][0], lattice[i][1], lattice[i][2]);

  /* print Element types and number of each element */
  fprintf(fp,"      %s         %s\n",elemM,elemX);
  fprintf(fp,"      %d         %d\n",n,2*n);
  fprintf(fp,"Cartesian\n");

  /* print coordinates for metal ions followed by coordinates for the 
     chalcogen ions */
  for (i = 0; i < n; i++)   fprintf(fp,"        %.9f                %.9f                %.9f\n",\
				    loc_m[i].x, loc_m[i].y, loc_m[i].z);
  for (i = 0; i < 2*n; i++) fprintf(fp,"        %.9f                %.9f                %.9f\n",\
				    loc_x[i].x, loc_x[i].y, loc_x[i].z);

  fclose(fp);
}

void print_help()
{
  print_version();

 fprintf(stdout,"Used to generate various sized unit cells for\n");
 fprintf(stdout,"1T and 2H trilayers of transition metal \n");
 fprintf(stdout,"dichalcogenides (TMD) with the MX2 structure.\n");
 fprintf(stdout,"\n");

 fprintf(stdout,"All options below must be specified in the\n");
 fprintf(stdout,"following order:\n");
 fprintf(stdout,"\t(1)  Lattice parameter a\n");
 fprintf(stdout,"\t(2)  Lattice parameter c\n");
 fprintf(stdout,"\t(3)  Super cell length a\'\n");
 fprintf(stdout,"\t(4)  Super cell length a\'\'\n");
 fprintf(stdout,"\t(5)  Super cell length b\'\n");
 fprintf(stdout,"\t(6)  Super cell length b\'\'\n");
 fprintf(stdout,"\t(7)  Layers (0 for Bulk)\n");
 fprintf(stdout,"\t(8)  1T (T) or 2H (F)\n");
 fprintf(stdout,"\t(9)  Randomize coordinates (T/F)\n");
 fprintf(stdout,"\t(10) Element M (use atomic number)\n");
 fprintf(stdout,"\t(11) Element X (use atomic number)\n");

 fprintf(stdout,"\n");
 fprintf(stdout,"The following options are optional, however\n");
 fprintf(stdout,"12-15 must be specified together and 16-17\n");
 fprintf(stdout,"must be specified together:\n");
 fprintf(stdout,"\t(12) Logical for Strain (T/F)\n");
 fprintf(stdout,"\t(13) Strain a axis (T/F)\n");
 fprintf(stdout,"\t(14) Strain b axis (T/F)\n");
 fprintf(stdout,"\t(15) Strain c axis (T/F)\n");
 fprintf(stdout,"\t(16) Minimum strain (%%)\n");
 fprintf(stdout,"\t(17) Maximum strain (%%)\n");
 fprintf(stdout,"\n");
}

void print_version()
{
 fprintf(stdout,"TMD CDW Unit Cell Generator %s\n",VERSION);
 fprintf(stdout,"Copyright (C) 2016-2017 David C. Miller\n");
 fprintf(stdout,"\n");
 fprintf(stdout,"Code written by David C. Miller (mill2723 at msu dot edu)\n");
 fprintf(stdout,"\n");
 fprintf(stdout,"This code is licensed under the GNU LGPL License\n");
 fprintf(stdout,"For details see <http://www.gnu.org/licenses/>\n");
 fprintf(stdout,"There is NO warranty; not even for MERCHANTABILITY or\n");
 fprintf(stdout,"FITNESS FOR A PARTICULAR PURPOSE.\n");
 fprintf(stdout,"\n");
}

void print_test_start()
{
  print_version();

  /*fprintf(stdout,"\n\n"); */
 fprintf(stdout,"---------------------------------------------------\n");
 fprintf(stdout,"Test Suite Output for TMD CDW Unit Cell Generator\n");
 fprintf(stdout,"---------------------------------------------------\n");
 fprintf(stdout,"Test Output Starting...\n");
 fprintf(stdout,"\n");
 fprintf(stdout,"Starting VASP Monolayer Test\n");
  test_print_vasp_monolayer();
 fprintf(stdout,"VASP Monolayer Test Output Complete (Comparison Pending...)\n");

 fprintf(stdout,"Starting VASP Bulk Test\n");
  test_print_vasp_bulk();
 fprintf(stdout,"VASP Bulk Test Output Complete (Comparison Pending...)\n");


}

/*
  deprecated function so precompiler statements to ignore
  comment or remove these if this function is needed
*/

#if 0 

void print_vasp\
(Location locTa[], Location locS[], unsigned n,\
 double lattice[3][3], char* name, char* elemM, char* elemX)
{
  unsigned i;
  fprintf(stdout,"%s\n",name);
  fprintf(stdout,"%.1f\n",1.0);

  /* the following output is formatted with spaces rather than tabs */
  for (i = 0; i < 3; ++i)
    fprintf(stdout,"   %.9f                %.9f                %.9f\n",	\
	    lattice[i][0], lattice[i][1], lattice[i][2]);

  fprintf(stdout,"     %s         %s\n",elemM,elemX);
  fprintf(stdout,"      %d         %d\n",n,2*n);
  fprintf(stdout,"Cartesian\n");

  for (i = 0; i < n; i++) fprintf(stdout,"        %.9f                %.9f                %.9f\n", \
				  locTa[i].x,locTa[i].y,locTa[i].z);

  for (i = 0; i < 2*n; i++) fprintf(stdout,"        %.9f                %.9f                %.9f\n", \
				    locS[i].x,locS[i].y,locS[i].z);
}

#endif 
