evaluate_it & orbit_it
======================
Homework 5, S. Shockley
-----------------------
***

Contents
========
|- data/
|  |- *.png
|  |- *.csv
|- homework6.pdf
|- Makefile
|- README
|- src/
|  |- lib/
|  |  |- *.h
|  |  |- leap.c
|  |  |- nbody_derivs.c
|  |  |- nbody_generator.c
|  |  |- nbody_print.c
|  |  |- nbody_util.c
|  |  |- *.c
|  |- nbody_it.c
|  |- plots.py

File Description
================
data/
-----
  This folder contains outputs from this homework. The .csv files are given in
  the format requested by the homework for single-line (two particles per line)
  output of the 2-body solution. Each is named including the mode used to gener-
  ate the solution and the eccentricity, and is created by the Makefile target
  `plots`. It also contains the plots generated by the python script and embeded
  in the homework PDF.

homework5.pdf
-------------
  This file is simply the answers to the homework questions with embedded pic-
  tures.

Makefile
--------
  Contains the code for the compilation of the executable and generation of
  plots.

README
------
  README-ception

src/
----
  This folder contains all of the source code for the project, listed below.
   - `nbody_it.c`: this is the code for the main and the loops for the two
     integration modes.
   - `plots.py`: this file turns the csv's inside the data/ folder into the plots
     for the homework.
  Then, in the lib subdirectory:
   - `*.h`: various header files. The int_methods.h and nbody.h files were writ-
     ten by myself, while the other two are NRIC's header files.
   - `leap.c`: the function that solves the leapfrog integration problem.
   - `nbody_derivs.c`: contains the two derivative functions for the leapfrog and
     RK4 integration methods.
   - `nbody_generator.c`: contains generation methods for the project, namely the
     2-body setup. Currently unimplemented are the random distributions.
   - `nbody_print.c`: intended for printing functions, currently unused.
   - `nbody_util.c`: contains utility functions such as vector indexing and
     reading functions.
   - `*.c`: all other files are NRIC functions.

Executables & make
==================

Makefile
--------
  Contains the following options:
  - `all`: compiles the executable.
  - `nbody_it`: eponymous.
  - `clean`: removes the executable.
  - `refresh`: removes and compiles the executable.
  - `plots`: generates the problem's plots.

nbody_it
-----------
  Solves an N-Body system. Input can be a generation method, and output can be
  to the console, a single file, or a directory printing multiple files. All
  outputs are given in the format `m x y z x* y* z*`. Commandline arguments of
  interest are:
  - `a`: the softening parameter. Defaults to 0.0.
  - `b`: two body mode. Default.
  - `e`: eccentricity for the two body mode. Must be between 0.0 and 1.0 (incl-
    usive).
  - `f`: file input mode. Argument must be the file name.
  - `h`: stepsize value. Must be positive.
  - `m`: integration mode. Must be "rk4" or "leap".
  - `n`: output frequency in steps per print. Must be positive and non-zero.
  - `o`: output directory mode. Prints a file per output in a given directory.
  - `O`: output file mode. Prints a line per step in a given filename.
  - `s`: steps. Must be positive.