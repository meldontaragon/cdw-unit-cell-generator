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
   tmd_functions.cpp
   code written by David Miller
   mill2723 at msu dot edu

   Additional functions to generator the coordinates for TMD trilayers.

***********************************************************************/

#include "header.h"

using namespace std;

void makeMSite\
(vector<Location>& atomsM, unsigned num, vector<double> & orig_lattice,\
 vector< vector<int> > & supercell, bool randomize, bool inversion,\
 bool monolayer)
{
  if (!(inversion) && !(monolayer))
    atomsM = vector<Location>(num*2,Location());
  else
    atomsM = vector<Location>(num,Location());
  vector< vector<double> > frac_loc = vector< vector<double> >(num,vector<double>(3,0.0));
  generateFracCoord(frac_loc, num, supercell);

  //cout << "DEBUG: Checking values..." << endl;
  //cout << orig_lattice[2] << endl;
  //cout << frac_loc[num-1][1] << endl;

  //cout << "DEBUG: Making actual coordinates..." << endl;
  double r1 = 0, r2 = 0;
  double scale = orig_lattice[0] * 0.03;

  srand(static_cast<unsigned>(time(NULL)));
  for (unsigned n = 0; n < num; n++)
    {
      //cout << "DEBUG: n=" << n << endl;
      if (randomize)
	{
	  r1 = (((1.0*rand()/RAND_MAX)*2.0)-1.0)*scale;
	  r2 = (((1.0*rand()/RAND_MAX)*2.0)-1.0)*scale;
	}
      //cout << "DEBUG: r1=" << r1 << " r2=" << r2 << endl;

      //location is n1*a1 + n2*a2 but a2=|a2|*0.5*x -sqrt(3)*0.5*|a2|

      if (!(inversion) && !(monolayer))
	{
	  atomsM[n] =
	    Location( (frac_loc[n][0] * orig_lattice[0])
		      - (frac_loc[n][1]*0.5*orig_lattice[0]) + r1,
		      (sqrt(3)*0.5*frac_loc[n][1] * orig_lattice[1]) +r2,
		      0.25 * orig_lattice[2]);
	  atomsM[n+num] =
	    Location( (frac_loc[n][0] * orig_lattice[0])
		      - (frac_loc[n][1]*0.5*orig_lattice[0]) + r1,
		      (sqrt(3)*0.5*frac_loc[n][1] * orig_lattice[1]) +r2,
		      0.75 * orig_lattice[2]);
	}
      else
	{
	  atomsM[n] =
	    Location( (frac_loc[n][0] * orig_lattice[0])
		      - (frac_loc[n][1]*0.5*orig_lattice[0]) + r1,
		      (sqrt(3)*0.5*frac_loc[n][1] * orig_lattice[1]) +r2,
		      0.5 * orig_lattice[2]);
	}
    }
}//end of makeMSite(...)

void makeXSite(vector<Location>& atomsX, unsigned num, vector<double> & orig_lattice,
	       vector< vector<int> > & supercell, bool inversion, bool monolayer)
{
  if (!(inversion) && !(monolayer))
    atomsX = vector<Location>(num*4,Location());
  else
    atomsX = vector<Location>(num*2,Location());

  //fractional location of atom within supercell
  vector< vector<double> > frac_loc = vector< vector<double> >(2*num,vector<double>(3,0.0));
  generateFracCoord(frac_loc, num, supercell);

  //correction vector for displacement of X atoms relative to M sites
  vector< vector<double> > correction = vector< vector<double> >(2,vector<double>(2,0.0));
  for (unsigned n = 0; n < num; n++)
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
      //atomsM[n] = Location((frac_loc[n][0] * orig_lattice[0]) - (frac_loc[n][1]*0.5*orig_lattice[0]) + r1,(sqrt(3)*0.5*frac_loc[n][1] * orig_lattice[1]) +r2, 0.5 * orig_lattice[2]);
      
      
      if (!(inversion) && !(monolayer))
	{
	  atomsX[n] =
	    Location( ((frac_loc[n][0] + correction[0][0]) * orig_lattice[0])
		      - ((frac_loc[n][1] + correction[0][1])*0.5*orig_lattice[0]),
		      (frac_loc[n][1] + correction[0][1]) * sqrt(3)*0.5*orig_lattice[1],
		      0.125 * orig_lattice[2]);
      
	  atomsX[n+num] =
	    Location( ((frac_loc[n][0] + correction[1][0]) * orig_lattice[0])
		      - ((frac_loc[n][1] + correction[1][1])*0.5*orig_lattice[0]),
		      (frac_loc[n][1] + correction[1][1]) * sqrt(3)*0.5*orig_lattice[1],
		      0.375 * orig_lattice[2]);
	  
	  atomsX[n+2*num] =
	    Location( ((frac_loc[n][0] - correction[0][0]) * orig_lattice[0])
		      - ((frac_loc[n][1] - correction[0][1])*0.5*orig_lattice[0]),
		      (frac_loc[n][1] - correction[0][1])
		      * sqrt(3)*0.5*orig_lattice[1],
		      0.625 * orig_lattice[2]);
	    
	  atomsX[n+3*num] =
	    Location( ((frac_loc[n][0] - correction[1][0]) * orig_lattice[0])
		      - ((frac_loc[n][1] - correction[1][1])*0.5*orig_lattice[0]),
		      (frac_loc[n][1] - correction[1][1]) * sqrt(3)*0.5*orig_lattice[1],
		      0.875 * orig_lattice[2]);
	}
      else
	{
	  atomsX[n] =
	    Location( ((frac_loc[n][0] + correction[0][0]) * orig_lattice[0])
		      - ((frac_loc[n][1] + correction[0][1])*0.5*orig_lattice[0]),
		      (frac_loc[n][1] + correction[0][1]) * sqrt(3)*0.5*orig_lattice[1],
		      0.75 * orig_lattice[2]);
      
	  atomsX[n+num] =
	    Location( ((frac_loc[n][0] + correction[1][0]) * orig_lattice[0])
		      - ((frac_loc[n][1] + correction[1][1])*0.5*orig_lattice[0]),
		      (frac_loc[n][1] + correction[1][1]) * sqrt(3)*0.5*orig_lattice[1],
		      0.25 * orig_lattice[2]);
	}
            
    }
}//end of makeXSite(...)

void printXYZ(vector<Location> & locTa, vector<Location>& locS, unsigned n)
{
  cout << 3*n << endl;
  cout << "c-cdw-TaS2 super cell SQRT(13)xSQRT(13)x1 (in angstroms)" <<endl;
  
  for (unsigned j = 0; j < n; j++)
    {
      printf("Ta        %.4f        %.4f        %.4f\n",
	     locTa[j][0], locTa[j][1],locTa[j][2]);
      printf("S        %.4f        %.4f        %.4f\n",
	     locS[j][0], locS[j][1],locS[j][2]);
      printf("S        %.4f        %.4f        %.4f\n",
	     locS[j+n][0], locS[j+n][1],locS[j+n][2]);
    }
}

void printVASP\
(vector<Location> & locTa, vector<Location> & locS, unsigned n,		\
 vector< vector<double> > & lattice, string name, string elemM = "Ta", string elemX = "S")
{
  cout << name << endl;
  cout << "1.0" << endl;

  cout << setprecision(10);

  for (unsigned i = 0; i < 3; i++)
    printf("   %.9f                %.9f                %.9f\n",
	   lattice[i][0], lattice[i][1], lattice[i][2]);
  
  cout << "     " << elemM << "         " << elemX << "\n";
  cout << "     " << n << "         " << 2*n << endl;
  cout << "Cartesian\n";
  
  for (unsigned i = 0; i < n; i++) printf("        %.9f                %.9f                %.9f\n",
				     locTa[i][0],locTa[i][1],locTa[i][2]);
  
  for (unsigned i = 0; i < 2*n; i++) printf("        %.9f                %.9f                %.9f\n",
				       locS[i][0],locS[i][1],locS[i][2]);
}

double dtor(double deg)
{
  return deg*PI/180.0;
}

double getLatticeVectorAngle(int n, int m)
{
  return acos((n - 0.5*m)/(sqrt(n*n + m*m - n*m*1.0))); //angle = arccos(a*b)/|a|*|b|
}

void generateFracCoord(vector< vector<double> > &frac_loc, unsigned num, vector< vector<int> > &supercell)
{
  unsigned xmax = static_cast<unsigned>(abs(supercell[0][0]) + abs(supercell[1][0]));
  unsigned ymax = static_cast<unsigned>(abs(supercell[0][1]) + abs(supercell[1][1]));

  unsigned count = 0;
  double min_angle = getLatticeVectorAngle(supercell[0][0], supercell[0][1]);
  //cout << "DEBUG:  num=" << num  << endl;
  //cout << "DEBUG: xmax=" << xmax << endl;
  //cout << "DEBUG: ymax=" << ymax << endl;
  //cout << "DEBUG: angle=" << min_angle << endl;
  //define vectors for the a and b lattices
  for (unsigned ii = 0; ii < xmax; ii++)
    {
      //cout << "\tDEBUG: ii=" << ii << endl;
      for (unsigned jj = 0; jj < ymax; jj++)
	{
	  //cout << "\t\tDEBUG: jj=" << jj << endl;
	  if (count >= num) { break; }	  
	  
	  if ( (ii == 0) && (jj == 0) )
	    {
	      for (unsigned kk = 0; kk < 3; kk++)
		{
		  frac_loc[0][kk] = 0;
		}
	      //cout << "\t\tCOORD[" << count << "]: (" << ii << ", " << jj << ")" << endl;
	      count++;
	    }
	  else
	    {
	      //cout << "\t\tCurrent Angle:" << getLatticeVectorAngle(ii,jj) << endl;
	      //cout << "\t\tAngle Param: " << ( getLatticeVectorAngle(ii,jj) - min_angle ) << endl;
	      if( ( getLatticeVectorAngle(static_cast<int>(ii),static_cast<int>(jj)) - min_angle ) > -EPS )
		{
		  //cout << "\t\tCOORD[" << count << "]: (" << ii << ", " << jj << ")" << endl;
		  frac_loc[count][0] = ii; ///xmax;
		  frac_loc[count][1] = jj; ///ymax;
		  frac_loc[count][2] = 0;
  		  count++;
		}
	    }	  
	}
    }
  if (count < num)
    {
      cerr << "Not enough coordinates made..." << endl;
      return;
    }
  //cout << "\"Fractional\" Coordinates Complete." << endl;
}

void printHelp()
{
  cout << "TMD CDW Unit Cell Generator"
       << endl
       << "Code written by David Miller "
       << "(mill2723 at msu dot edu)"
       << endl
       << "This code is licensed under the GNU LGPL License"
       << endl
       << "For details see <http://www.gnu.org/licenses/>"
       << endl
       << "There is NO warranty; not even for MERCHANTABILITY or "
       << "FITNESS FOR A PARTICULAR PURPOSE."
       << endl
       << endl;
  cout << "Used to generate various sized unit cells for\n"
       << "1T and 2H trilayers of transition metal \n"
       << "dichalcogenides (TMD) with the MX2 structure."
       << endl
       << endl
       << endl;
  cout << "All options below must be specified in the following order:"
       << endl;
  cout << "\t(1)  Lattice parameter a"
       << endl;
  cout << "\t(2)  Lattice parameter c"
       << endl;
  cout << "\t(3)  Super cell length a\'"
       << endl;
  cout << "\t(4)  Super cell length a\'\'"
       << endl;
  cout << "\t(5)  Super cell length b\'"
       << endl;
  cout << "\t(6)  Super cell length b\'\'"
       << endl;
  cout << "\t(7)  Monolayer (T) or Bulk (F)"
       << endl;
  cout << "\t(8)  1T (T) or 2H (F)"
       << endl;
  cout << "\t(9)  Randomize coordinates (T/F)"
       << endl;
  cout << "\t(10) Element M (use atomic number)"
       << endl;
  cout << "\t(11) Element X (use atomic number)"
       << endl;  
}
