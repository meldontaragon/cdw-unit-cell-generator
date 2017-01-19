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
#define VERSION "0.0.0"
#endif
/* VERSION DO NOT EDIT */

#include "header.h"
void print_xyz(Location loc_m[], Location loc_x[], unsigned n)
{
  unsigned j;
  printf("%d\n",3*n);
  printf("c-cdw-TaS2 super cell SQRT(13)xSQRT(13)x1 (in angstroms)\n");

  for (j = 0; j < n; j++)
    {
      printf("Ta        %.4f        %.4f        %.4f\n",\
	     loc_m[j].x, loc_m[j].y,loc_m[j].z);
      printf("S        %.4f        %.4f        %.4f\n",\
	     loc_x[j].x, loc_x[j].y,loc_x[j].z);
      printf("S        %.4f        %.4f        %.4f\n",\
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
  printf("TMD CDW Unit Cell Generator %s\n",VERSION);
  printf("Copyright (C) 2016-2017 David C. Miller\n");
  printf("Code written by David C. Miller (mill2723 at msu dot edu)\n");
  printf("\n");
  printf("This code is licensed under the GNU LGPL License\n");
  printf("For details see <http://www.gnu.org/licenses/>\n");
  printf("There is NO warranty; not even for MERCHANTABILITY or\n");
  printf("FITNESS FOR A PARTICULAR PURPOSE.\n");
  printf("\n");

  printf("Used to generate various sized unit cells for\n");
  printf("1T and 2H trilayers of transition metal \n");
  printf("dichalcogenides (TMD) with the MX2 structure.\n");
  printf("\n");

  printf("All options below must be specified in the\n");
  printf("following order:\n");
  printf("\t(1)  Lattice parameter a\n");
  printf("\t(2)  Lattice parameter c\n");
  printf("\t(3)  Super cell length a\'\n");
  printf("\t(4)  Super cell length a\'\'\n");
  printf("\t(5)  Super cell length b\'\n");
  printf("\t(6)  Super cell length b\'\'\n");
  printf("\t(7)  Layers (0 for Bulk)\n");
  printf("\t(8)  1T (T) or 2H (F)\n");
  printf("\t(9)  Randomize coordinates (T/F)\n");
  printf("\t(10) Element M (use atomic number)\n");
  printf("\t(11) Element X (use atomic number)\n");

  printf("\n");
  printf("The following options are optional, however\n");
  printf("12-15 must be specified together and 16-17\n");
  printf("must be specified together:\n");
  printf("\t(12) Logical for Strain (T/F)\n");
  printf("\t(13) Strain a axis (T/F)\n");
  printf("\t(14) Strain b axis (T/F)\n");
  printf("\t(15) Strain c axis (T/F)\n");
  printf("\t(16) Minimum strain (%%)\n");
  printf("\t(17) Maximum strain (%%)\n");
  printf("\n");
}

/*
  deprecated function so precompiler statements to ignore
  comment or remove these if this function is needed
*/

#if 0

void printVASP\
(Location locTa[], Location locS[], unsigned n,\
 double lattice[3][3], char* name, char* elemM, char* elemX)
{
  unsigned i;
  printf("%s\n",name);
  printf("%.1f\n",1.0);

  /* the following output is formatted with spaces rather than tabs */
  for (i = 0; i < 3; ++i)
    printf("   %.9f                %.9f                %.9f\n",\
	   lattice[i][0], lattice[i][1], lattice[i][2]);

  printf("     %s         %s\n",elemM,elemX);
  printf("      %d         %d\n",n,2*n);
  printf("Cartesian\n");

  for (i = 0; i < n; i++) printf("        %.9f                %.9f                %.9f\n",\
				 locTa[i].x,locTa[i].y,locTa[i].z);

  for (i = 0; i < 2*n; i++) printf("        %.9f                %.9f                %.9f\n",\
				   locS[i].x,locS[i].y,locS[i].z);
}

#endif
