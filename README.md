# CDW Unit Cell Generator

Code written by David C. Miller (`mill2723 at msu dot edu`)

Michigan State University

Department of Physics and Astronomy

## Quick Install:
 
Try the following:
+ change INSTALLDIR (typically a folder in PATH like
  /usr/local/bin)
  - this is defined in makefile.include or can be passed
    through make
+ make install
  - this will run `make reset`, `make all`, and then
  install the binary files in the INSTALLDIR

Other options:
+ make dir
  - run this or manually make the necessary directories (`./bin`
  and `./obj`)
+ make final
  - final binary files will be in `./bin`
+ make link
  - will create symbolic links in INSTALLDIR for the binary files
+ make all
  - will make binary files using the debug optimization options
  (for distributed versions there is a chance that `make all`
  runs `make final` and `make debug` runs this case)
+ make uninstall
  - removes binary files from INSTALLDIR (whether linked or copied)

You may need to run any of the above as root depending on 
local folder privileges.

## Information

### Description
This project is developed to generate a VASP coordinate file for a
transition metal dichalcogenide (TMD) structure. Given a choice of
metal and chalcogen atoms, the lattice parameters, and the super cell
size in additional to information regarding the phase and layers of
the material.

### Current Functionality
Currently the only options for layers are bulk and monolayer systems.
In the future I hope to add functionallity to have an arbitrary number
of layer (i.e. 0 = bulk, 1 = monolayer, 2 = bilayer, etc.).

### Use
Use of this program involves calling the binary and giving it all the
arguments necessary for output. The order and style of inputs is given
in the help for the program (just call the binary without any arguments).
They are also given below:

All options below must be specified in the following order:

	1.  Lattice parameter a (in Angstroms)
	2.  Lattice parameter c (in Angstroms)
	3.  Super cell length a'
	4.  Super cell length a''
	5.  Super cell length b'
	6.  Super cell length b''
	7.  Monolayer (T) or Bulk (F)
	8.  1T (T) or 2H (F)
	9.  Randomize coordinates (T/F)
	10. Element M (use atomic number)
	11. Element X (use atomic number)

The a', a'', b', and b'' parameters are the two-vectors that specify the
super-cell size. These describe how the super-cell is constructed in terms
of the original a and b lattice parameters (where b is at 120 degrees to
a for the hexagonal primitive cell). This allows for not only typical super
cells. For example a 3x3 super cell would be (3,0) x (0,3). More complex
cells like the Star of David distortion in TaS2 are also possible where the
Sqrt(13)xSqrt(13) cell is given by (4,1) x (-1,3) which results in the proper
angle between the super cell and the primitive cell.

### Bugs

Please report bugs to me (`mill2723 at msu dot edu`) with the subject
"TMD CDW Coordinate Generator - Bug" and a description of the relevant
information about how to replicate the bug. Please include: the program
being run, version information, an attached input and output file (if
possible), the results of stdout and stderr, and a description of what
you expected. 

I will try to take a look at it as soon as possible and get back to
you with a fix or a patch.

### Compiling and Installing

This program was designed to be compiled with the GNU C Compiler
(gcc) and GNU Make. It also uses some simple bash commands. It *should*
work on any Linux/UNIX system that can meet these requirements but I can
offer no guarantees. I suspect (but have not tested) that this could be
compiled on Cygwin or any other system with a C compiler if one so
desired. If another C compiler is used it is very likely that the current
compiler options will need to be changed. This program does not use
any special libraries so it should suffice to pass the include and source
directories (`./src` and `./inc`).

Currently I have used this on Ubuntu 16.04, and CentOS 6.6.

Dependencies:

+ GNU gcc 4.4.7
+ GNU Make 3.81 (again, this is what I have run it with but I suspect
it would work with an earlier version)
+ GNU bash 4.1.2 (same as above)

Also, this does compile with the C90 ISO standards. This includes
using the gcc `-ansi` option. If you wish to compile without this
dependency there should be no problems. This is mostly to ensure
better portability. Additionally, one may wish to remove the `WARN_FLAG`
variable from the makefile or at least turnoff `-Werror`.

Also, the code was written on an x64 system but there shouldn't be
many dependencies if any. If you do compile and run this code 
on a widely different system (or perhaps just a very popular system)
than what has been described please let me know so I can add it to
this document.

To install:

(a) Modify the makefile.include file to match your desired
settings. This may include changing the compiler, flags,
install directory, or programs to be compiled. Again,
I can make no guareentees about the operation of the
program without using the default settings on a system
that it has been tested on.

(b) Run `make`

(c) Assuming everything went smoothly run `make install`
(you may need to run this as root or use sudo depen-
ding on the installation directory). I was also ver-
ify that there is nothing in the install directory
with an identical name to one of the binaries being
added. The makefile *will* overwrite it.

-  To uninstall simply run `make uninstall`
-  The install script only makes a link from the binary
directory to the bin/ directory for the local compi-
lation. Thus, if you recompile, there is no need to
reinstall. If you would rather the binary files are
copied, use `make final`

Other: `make clean` will delete any extraneous files.
`make reset` will delete the binary files.

### Git repository

This code is available at:
https://github.com/david-c-miller/cdw-unit-cell-generator.git

### Disclaimer

I wrote this to be a useful program for computational work 
done with VASP. There are other tools out there as well and
I cannot promise compatability with them or that my code will
produce identical results. Please verify the results and 
output on your own.

### License

This code is licensed under LGPLv3. See LICENSE for details.
