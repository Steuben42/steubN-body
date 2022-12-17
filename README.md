nbody_it
======================
S. Shockley
-----------------------
***

TO RUN THIS ON YOUR OWN, IT IS RECOMMENDED TO USE THE DEFAULT OPTS FILE 
PROVIDED IN THE OPTS FOLDER, AND TO USE THE SOLAR GENERATOR FILE. YOU
MUST FIRST MAKE THE EXECUTABLES, BY SIMPLY USING `make all`, AND THEN
GENERATE THE INITIAL CONDITIONS VIA 
`./generate_nbody -f ./opts/objs/solar -o ./data/{init_file.csv}`
THEN RUNNING THE PROGRAM BY COPYING THE DEFAULT SETTINGS INTO A NEW 
FILE `{opts_file}` AND SETTING `file_in:{opts_file}`, THEN RUNNING
`./nbody_it -F ./opts/{opts_file}`.

Contents
========
|- data/
|  |- .gitignore
|- docs/
|  |- imgs/
|  |  |- *
|  |- Shockley_Term_Project.pdf
|  |- term_project_ASTR415.pptx
|- Makefile
|- opts/
|  |- default
|  |- objs/
|  |  |- solar
|- README.md
|- src/
|  |- generate_nbody.c
|  |- lib/
|  |  |- *.h
|  |  |- bhmode_nbody.c
|  |  |- derivs_nbody.c
|  |  |- generate_nbody.c
|  |  |- leap.c
|  |  |- print_nbody.c
|  |  |- util_.nbody.c
|  |  |- *.c
|  |- nbody_it.c
|  |- python/
|  |  |- *.py
|  |  |- trojans_plot.py

File Description
================
data/
-----
  This folder is intended to be used for input and output data of the 
  simulation, i.e., all of the generated and inputted space-delimited .csv
  files.

docs/
-------------
  This file contains the `imgs/` folder, in which are all the generated plots
  for this project, most of which are included in the first file, the .pdf.
  The .pdf provides an overview of the code and the assignment and is the
  writeup for the term project. The next file, the .pptx, are the slides used
  in the Dec. 8th presentation of the code.

Makefile
--------
  Contains the code for the compilation of the executables.

opts/
-----
  This file contains the automatically generated `default` script settings file,
  and a subdirectory `objs/` which contains an example object script. To recreate
  the simulation, refer to the .pdf from the `docs/` folder for options. Namely:
  do not use the `@GEN` scope, and use the `file_in:` parameter using a generated
  file from the generator executable which is fed files such as `objs/`. Although
  there is some code for a third structure, `@SATELLITE`, results are unreliable,
  and so it is recommended to exclusively use the `@GROUP` and `@RING` 
  declarations.

README
------
  README-ception

src/
----
  This folder contains all of the source code for the project, listed below.
   - `generate_nbody.c`: this is the code for a helpful initial conditions
     generator tool which is used for processing the objs/ files.
   - `nbody_it.c`: this is the code for the main and original loop.
  Then, in the lib subdirectory:
   - `*.h`: various header files. The int_methods.h and nbody.h files were writ-
     ten by myself, while the other two are NRIC's header files.
   - `bhmode_nbody.c`: the code for the implementation of the Barnes-Hut 
     algorithm, included tree building and manipulation and structure
     intialization.
   - `derivs_nbody.c`: contains the derivs functions for the non-BH N-body method
     and the BH implementation.
   - `generator_nbody.c`: contains generation methods for the project.
   - `leap.c`: the function that solves the leapfrog integration problem.
   - `print_nbody.c`: printing functions primarily used by the generator.
   - `util_nbody.c`: contains utility functions such as vector indexing and
     reading functions as well as type conversion.
   - `*.c`: all other files are NRIC functions.
  Lastly, in the `python/` subdirectory:
   - `trojans_plot.py`: the code used to make the images of the Jupiter-belt
     simulation. To adapt, change the `path` string to the subdirectory
     the files are contained in and the `range()` value to the number of frames.
     Note: this assumes that the output frequency is 10, and that BH-mode was
     used.

Executables & make
==================

Makefile
--------
  Contains the following options:
  - `all`: compiles the executables.
  - `nbody_it`: eponymous. Generates default file.
  - `generate_nbody`: eponymous.
  - `clean`: removes the executables.
  - `refresh`: removes and compiles the executables.
  - `backup`: intended for automatic backups.

generate_nbody
--------------
  Generates an initial conditions file for the program. Two commandline 
  arguments will likely be used:
  - `-f`: the input object script.
  - `-o`: the output space delimited (.csv) file.

nbody_it
-----------
  Solves an N-Body system. Input can be a generation method, and output can be
  to the console, a single file, or a directory printing multiple files. All
  outputs are given in the format `m x y z x* y* z*`. The commandline argument
  of interest is:
  - `F`: use an input options script.
  In addition, one can pass `-g {N}` while using the sphere generation mode,
  which will specify the number of puts to generate in the sphere.