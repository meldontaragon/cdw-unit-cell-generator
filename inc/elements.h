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
/* **********************************************
   elements.h
   code written by David Miller
   mill2723 at msu dot edu

   Contains definitions for lists of all the elements 
   in the periodic table. 
**************************************************/

#ifndef ELEMENTS_H
#define ELEMENTS_H

static const unsigned MAX_ELEMENT = 117;

static const char * Symbol[] = {"H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Uut","Fl","Uup","Lv","Uus","Uuo"};

/*
static const char * Name[] = {"Hydrogen","Helium","Lithium","Beryllium","Boron","Carbon","Nitrogen","Oxygen","Fluorine","Neon","Sodium","Magnesium","Aluminum","Silicon","Phosphorus","Sulfur","Chlorine","Argon","Potassium","Calcium","Scandium","Titanium","Vanadium","Chromium","Manganese","Iron","Cobalt","Nickel","Copper","Zinc","Gallium","Germanium","Arsenic","Selenium","Bromine","Krypton","Rubidium","Strontium","Yttrium","Zirconium","Niobium","Molybdenum","Technetium","Ruthenium","Rhodium","Palladium","Silver","Cadmium","Indium","Tin","Antimony","Tellurium","Iodine","Xenon","Cesium","Barium","Lanthanum","Cerium","Praseodymium","Neodymium","Promethium","Samarium","Europium","Gadolinium","Terbium","Dysprosium","Holmium","Erbium","Thulium","Ytterbium","Lutetium","Hafnium","Tantalum","Tungsten","Rhenium","Osmium","Iridium","Platinum","Gold","Mercury","Thallium","Lead","Bismuth","Polonium","Astatine","Radon","Francium","Radium","Actinium","Thorium","Protactinium","Uranium","Neptunium","Plutonium","Americium","Curium","Berkelium","Californium","Einsteinium","Fermium","Mendelevium","Nobelium","Lawrencium","Rutherfordium","Dubnium","Seaborgium","Bohrium","Hassium","Meitnerium","Darmstadtium","Roentgenium","Copernicium","Ununtrium","Flerovium","Ununpentium","Livermorium","Ununseptium","Ununoctium"};
*/
enum AtomicSymbol_E { H=1 , He , Li , Be , B , C , N , O , F , Ne , Na , Mg , Al , Si , P , S , Cl , Ar , K , Ca , Sc , Ti , V , Cr , Mn , Fe , Co , Ni , Cu , Zn , Ga , Ge , As , Se , Br , Kr , Rb , Sr , Y , Zr , Nb , Mo , Tc , Ru , Rh , Pd , Ag , Cd , In , Sn , Sb , Te , I , Xe , Cs , Ba , La , Ce , Pr , Nd , Pm , Sm , Eu , Gd , Tb , Dy , Ho , Er , Tm , Yb , Lu , Hf , Ta , W , Re , Os , Ir , Pt , Au , Hg , Tl , Pb , Bi , Po , At , Rn , Fr , Ra , Ac , Th , Pa , U , Np , Pu , Am , Cm , Bk , Cf , Es , Fm , Md , No , Lr , Rf , Db , Sg , Bh , Hs , Mt , Ds , Rg , Cn , Uut , Fl , Uup , Lv , Uus , Uuo };

enum  AtomicName_E { Hydrogen , Helium , Lithium , Beryllium , Boron , Carbon , Nitrogen , Oxygen , Fluorine , Neon , Sodium , Magnesium , Aluminum , Silicon , Phosphorus , Sulfur , Chlorine , Argon , Potassium , Calcium , Scandium , Titanium , Vanadium , Chromium , Manganese , Iron , Cobalt , Nickel , Copper , Zinc , Gallium , Germanium , Arsenic , Selenium , Bromine , Krypton , Rubidium , Strontium , Yttrium , Zirconium , Niobium , Molybdenum , Technetium , Ruthenium , Rhodium , Palladium , Silver , Cadmium , Indium , Tin , Antimony , Tellurium , Iodine , Xenon , Cesium , Barium , Lanthanum , Cerium , Praseodymium , Neodymium , Promethium , Samarium , Europium , Gadolinium , Terbium , Dysprosium , Holmium , Erbium , Thulium , Ytterbium , Lutetium , Hafnium , Tantalum , Tungsten , Rhenium , Osmium , Iridium , Platinum , Gold , Mercury , Thallium , Lead , Bismuth , Polonium , Astatine , Radon , Francium , Radium , Actinium , Thorium , Protactinium , Uranium , Neptunium , Plutonium , Americium , Curium , Berkelium , Californium , Einsteinium , Fermium , Mendelevium , Nobelium , Lawrencium , Rutherfordium , Dubnium , Seaborgium , Bohrium , Hassium , Meitnerium , Darmstadtium , Roentgenium , Copernicium , Ununtrium , Flerovium , Ununpentium , Livermorium , Ununseptium , Ununoctium };

typedef enum AtomicSymbol_E AtomicSymbol;
typedef enum AtomicName_E AtomicName;
#endif
