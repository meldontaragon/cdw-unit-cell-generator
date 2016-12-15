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
   tmd_arr_functions.c
   code written by David Miller
   mill2723 at msu dot edu

   Additional functions to generator the coordinates for TMD trilayers.
   All of these functions use arrays as opposed to vectors.

***********************************************************************/

#include "header_c.h"
#include <math.h>
#include <time.h>

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

void printXYZ(Location locTa[], Location locS[], unsigned n)
{
  unsigned j;
  printf("%d\n",3*n);
  printf("c-cdw-TaS2 super cell SQRT(13)xSQRT(13)x1 (in angstroms)\n");
  
  for (j = 0; j < n; j++)
    {
      printf("Ta        %.4f        %.4f        %.4f\n",\
	     locTa[j].x, locTa[j].y,locTa[j].z);
      printf("S        %.4f        %.4f        %.4f\n",\
	     locS[j].x, locS[j].y,locS[j].z);
      printf("S        %.4f        %.4f        %.4f\n",\
	     locS[j+n].x, locS[j+n].y,locS[j+n].z);
    }
}

void printVASP\
(Location locTa[], Location locS[], unsigned n,\
 double lattice[3][3], char* name, char* elemM, char* elemX)
{
  unsigned i;
  printf("%s\n",name);
  printf("%.1f\n",1.0);

  for (i = 0; i < 3; i++)
    printf("   %.9f                %.9f                %.9f\n",\
	   lattice[i][0], lattice[i][1], lattice[i][2]);

  printf("     %s         %s\n",elemM,elemX);
  printf("      %d         %d\n",n,2*n);
  printf("Cartesian\n");
  
  for (i = 0; i < n; i++) printf("        %.9f                %.9f                %.9f\n",\
				 locTa[i].x,locTa[i].y,locTa[i].z);
  
  for (i = 0; i < 2*n; i++) printf("        %.9f                %.9f                %.9f\n",\
				   locS[i].x,locS[i].y,locS[i].x);
}

double dtor(double deg)
{
  return deg*PI/180.0;
}

double getLatticeVectorAngle(int n, int m)
{
  /* angle = arccos(a*b)/|a|*|b| */
  return acos((n - 0.5*m)/(sqrt(n*n + m*m - n*m*1.0))); 
}

void generateFracCoord(double frac_loc[][3], unsigned num, int supercell[2][2])
{
  unsigned xmax, ymax, count, ii, jj, kk;
  double min_angle;
  
  xmax = (unsigned)(abs(supercell[0][0]) + abs(supercell[1][0]));
  ymax = (unsigned)(abs(supercell[0][1]) + abs(supercell[1][1]));

  count = 0;
  min_angle = getLatticeVectorAngle(supercell[0][0], supercell[0][1]);

  for (ii = 0; ii < xmax; ii++)
    {
      for (jj = 0; jj < ymax; jj++)
	{
	  if (count >= num) { break; }	  
	  
	  if ( (ii == 0) && (jj == 0) )
	    {
	      for (kk = 0; kk < 3; kk++)
		{
		  frac_loc[0][kk] = 0;
		}
	      count++;
	    }
	  else
	    {
	      if( ( getLatticeVectorAngle((int)(ii),(int)(jj)) - min_angle ) > -EPS )
		{
		  frac_loc[count][0] = ii; 
		  frac_loc[count][1] = jj; 
		  frac_loc[count][2] = 0;
  		  count++;
		}
	    }	  
	}
    }
  if (count < num)
    {
      printf("Not enough coordinates made...\n");
      return;
    }
}

void printHelp()
{
  printf("TMD CDW Unit Cell Generator");
  printf("\n");
  printf("Code written by David Miller ");
  printf("(mill2723 at msu dot edu)");
  printf("\n");
  printf("This code is licensed under the GNU LGPL License");
  printf("\n");
  printf("For details see <http://www.gnu.org/licenses/>");
  printf("\n");
  printf("There is NO warranty; not even for MERCHANTABILITY or ");
  printf("FITNESS FOR A PARTICULAR PURPOSE.");
  printf("\n");
  printf("\n");;
  printf("Used to generate various sized unit cells for\n");
  printf("1T and 2H trilayers of transition metal \n");
  printf("dichalcogenides (TMD) with the MX2 structure.");
  printf("\n");
  printf("\n");
  printf("\n");;
  printf("All options below must be specified in the following order:");
  printf("\n");;
  printf("\t(1)  Lattice parameter a");
  printf("\n");;
  printf("\t(2)  Lattice parameter c");
  printf("\n");;
  printf("\t(3)  Super cell length a\'");
  printf("\n");;
  printf("\t(4)  Super cell length a\'\'");
  printf("\n");;
  printf("\t(5)  Super cell length b\'");
  printf("\n");;
  printf("\t(6)  Super cell length b\'\'");
  printf("\n");;
  printf("\t(7)  Monolayer (T) or Bulk (F)");
  printf("\n");
  printf("\t(8)  1T (T) or 2H (F)");
  printf("\n");
  printf("\t(9)  Randomize coordinates (T/F)");
  printf("\n");
  printf("\t(10) Element M (use atomic number)");
  printf("\n");
  printf("\t(11) Element X (use atomic number)");
  printf("\n");  
}

int atob(char a)
{
  if ( (a == 't') || (a == 'T') ) return 1;
  else if ( (a == 'f') || (a == 'F') ) return 0;
  else exit(-1);
}
