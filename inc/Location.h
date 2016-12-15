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

#include <stdio.h>
  
void setCoord(Location *loc, const double ix, const double iy, const double iz)
{
  loc->x = ix; loc->y = iy; loc->z = iz;
}

void setElement(Location *loc, const AtomicSymbol ielem)
{
  if (ielem > MAX_ELEMENT) {
    printf("Element given is outside known range\n");
  }
  else
    {
      loc->elem = ielem;
    }
}

#endif
