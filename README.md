# CDW Unit Cell Generator

Copyright (C) 2016-2017 David C. Miller

Code written by David C. Miller (`mill2723 at msu dot edu`)

Michigan State University

Department of Physics and Astronomy

## Quick Install:
 
Try the following:

+ `./configure`
+ `make`
+ `make install`

Other options:

+ `make debug`
  - runs make with '-O0' optimization flag and '-pg -g3' debug flags
+ `make optim`
  - runs make with '-O3' optimization flag
+ `make uninstall`
  - removes binary files from install directory

You may need to run any of the above as root depending on local folder
privileges. Additonally, see `./configure --help` for more options.

## Information

### Description
This project is developed to generate a VASP coordinate file for a transition
metal dichalcogenide (TMD) structure. Given a choice of metal and chalcogen
atoms, the lattice parameters, and the super cell size in additional to
information regarding the phase and layers of the material.

### Current Functionality
Currently the only options for layers are bulk and monolayer systems. In the
future I hope to add functionallity to have an arbitrary number of layer (i.e. 0
for bulk, 1 for monolayer, 2 for bilayer, etc.). 

### Use
Use of this program involves calling the binary and giving it all the arguments
necessary for output. The order and style of inputs is given in the help for the
program (just call the binary without any arguments). They are also given below:

All of the options in the first section below must be specified in the following
order: 

       1. Lattice parameter a (in Angstroms)
       2. Lattice parameter c (in Angstroms)
       3. Super cell length a'
       4. Super cell length a''
       5. Super cell length b'
       6. Super cell length b''
       7. Layers (Bulk = 0, Monolayer = 1)
       8. 1T (T) or 2H (F)
       9. Randomize coordinates (T/F)
       10. Element M (use atomic number)
       11. Element X (use atomic number)

The following options are not required but allow generation of additional unit
cells with strain added along different axes:

      12. Logical for Strain (T/F)
      13. Strain a axis (T/F)
      14. Strain b axis (T/F)
      15. Strain c axis (T/F)

In addition to those options above, the following two options can be added after
to specify the amount of strain:
   
      16. Minimum strain (%)
      17. Maximum strain (%)

The a', a'', b', and b'' parameters are the two-vectors that specify the
super-cell size. These describe how the super-cell is constructed in terms of
the original a and b lattice parameters (where b is at 120 degrees to a for the
hexagonal primitive cell). This allows for not only typical super cells. For
example a 3x3 super cell would be (3,0) x (0,3). More complex cells like the
Star of David distortion in TaS2 are also possible where the $\sqrt{13}x\sqrt{13}$
cell is given by (4,1) x (-1,3) which results in the proper angle between the
super cell and the primitive cell.



### Compiling and Installing

This program was designed to be compiled with the GNU C Compiler (gcc) and GNU
Make. Additionally  It can also use the clang processor. It also uses some simple bash
commands. It *should* work on any Linux/UNIX system that can meet these
requirements but I can offer no guarantees. I suspect (but have not tested) that
this could be compiled on Cygwin or any other system with a C compiler if one so
desired. If another C compiler is used it is very likely that the current
compiler options will need to be changed. This program does not use any special
libraries so it should suffice to pass the include and source directories (./src
and ./inc). 

Currently I have used this on Ubuntu 16.04, and CentOS 6.6.

Dependencies (with earliest tested version):

+ GNU gcc 4.6.4
+ GNU Make 3.81
+ GNU bash 4.1.2 
+ GNU automake 1.11.1
+ GNU autoconf 2.63

Also, this does compile with the C99 ISO standards. This includes using the gcc
`-std=c99` option. If you wish to compile without this dependency there should be
no problems. This is mostly to ensure better portability. Additionally, one may
wish to remove the `WARN_FLAG` variable from the makefile.

Also, the code was written on an x64 system but there shouldn't be many
dependencies if any. If you do compile and run this code on a widely different
system (or perhaps just a very popular system) than what has been described
please let me know so I can add it to this document.

It should be possible to compile this code with older version of gcc or with
Clang by removing certain warning flags.

To install:

1. Run `./configure` with any additonal options (see `./configure --help` for
details). Use `--prefix=/path/to/install/root` to change the root directory for
the install. By default `make install` will install the binary file to
'/usr/local/bin'.
2. Run `make`
3. Assuming everything went smoothly run `make install` (you may need to run
this as root or use sudo depending on the installation directory). 

To uninstall simply run `make uninstall`

Other: `make mostlyclean` will delete any extraneous files while `make clean`
will remove all files generated by `make`

### Git repository

This code is available at:
[repository](https://github.com/david-c-miller/cdw-unit-cell-generator)

For releases see:
[releases](https://github.com/david-c-miller/cdw-unit-cell-generator/releases/latest)

For other information go to:
[cdw-unit-cell-generator](https://david-c-miller.github.io/cdw-unit-cell-generator/)

### Bug Reports

Please reports all bugs in the GitHub
[Issues](https://github.com/david-c-miller/cdw-unit-cell-generator/issues)
section. This way others can view the bugs. If there is no response there you
can try sending me an email with the subject "Bug - CDW Unit Cell Generator"
along with a thorough description of the bug. Please include the following:

+ output files
+ full printout of stdout and stderr
+ version
+ any other potentially useful information
+ description of what you expected and what went wrong

I will try to take a look at it as soon as possible and get back to you with a
fix or a patch.

### Disclaimer

I wrote this to be a useful program for computational work done with VASP. There
are other tools out there as well and I cannot promise compatability with them
or that my code will produce identical results. Please verify the results and
output on your own. Additionally, this project is partially done in my free time
so I make no guarantees about the time line for fixes or additional
features. Likely the best way to fix something is to modify the code yourself
and make a pull request.

### License

This code is licensed under LGPLv3. There is NO warranty; not even for
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See LICENSE for details.

