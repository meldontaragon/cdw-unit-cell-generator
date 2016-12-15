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
/* *****************************************************
   Location.h
   Code written by David Miller
   mill2723 at msu dot edu

   This contain the definition for a class that holds
   the location of an atom at a cartesian grid location
   with a given atomic symbol
*/

#ifndef LOCATION_H
#define LOCATION_H

#include "header.h"
#include "elements.h"

using namespace std;

class Location
{
 private:
  double x, y, z;
  Elements::AtomicSymbol elem;
  
 public:
  Location() 
    {
      x = 0;
      y = 0;
      z = 0;
      elem = Elements::H;
    }

  Location(const double ix, const double iy, const double iz)
    {
      x = ix;
      y = iy;
      z = iz;
      elem = Elements::H;
    }

  Location(const double ix, const double iy, const double iz, const Elements::AtomicSymbol ielem)
    {
      x = ix;
      y = iy;
      z = iz;
      if (ielem > Elements::MAX_ELEMENT) {
	cerr << "Element given is outside known range" << endl;
	elem = Elements::H;
      }
      else
	elem = ielem;
    }

  Location(const vector<double> &ic)
    {
      if (ic.size() != 3)
	{
	  cerr << "Size of coordinate vector is not 3." << endl;
	  x = 0; y = 0; z = 0;
	}
      else
	{
	  x = ic[0]; y = ic[1]; z = ic[2];
	}
      elem = Elements::H; 
    }

  Location(const vector<double> &ic, const Elements::AtomicSymbol ielem)
    {
      if (ic.size() != 3)
	{
	  cerr << "Size of coordinate vector is not 3." << endl;
	  x = 0; y = 0; z = 0;
	}
      else
	{
	  x = ic[0]; y = ic[1]; z = ic[2];
	}
      if (ielem > Elements::MAX_ELEMENT) {
	cerr << "Element given is outside known range" << endl;
	elem = Elements::H;
      }
      else
	elem = ielem;
    }

  double &operator[](const int n)
  {
    switch(n)
      {
      case 0:
	return x;
      case 1:
	return y;
      case 2:
	return z;
      default:
	cout << "Error: Index out of bounds\n";
	throw -1;
      }
  }
  
  void setCoord(const double ix, const double iy, const double iz)
  {
    x = ix; y = iy; z = iz;
  }
  void setCoord(const vector<double> &ic)
  {
    if (ic.size() != 3)
      {
	cerr << "Size of coordinate vector is not 3." << endl;
	// x = 0; y = 0; z = 0;
      }
    else
      {
	x = ic[0]; y = ic[1]; z = ic[2];
      }
  }

  vector<double> getCoord()
    {
      vector<double> tmp = vector<double>(3,0);
      tmp[0] = x;
      tmp[1] = y;
      tmp[2] = z;
      return tmp;
    }

  double getX()
  {
    return x;
  }
  double getY()
  {
    return y;
  }
  double getZ()
  {
    return z;
  }

  string getElementSymbol()
  {
    return Elements::Symbol[elem];
  }
  Elements::AtomicSymbol getElement()
    {
      return elem;
    }

  void setElement(const Elements::AtomicSymbol ielem)
  {
    if (ielem > Elements::MAX_ELEMENT) {
      cerr << "Element given is outside known range" << endl;
      //elem = 0;
    }
    else
      {
	elem = ielem;
      }
  }
};

#endif
