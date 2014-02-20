This is SuperNu, a IMC+DDMC Radiation Transport Code.


SETUP INSTRUCTIONS
==================
The following procedure should produce a working supernu executable.

- if you do not have the hg-versioninfo module installed:
  # cp ./Tools/hgversioninfo.py ~/py/
  # export PYTHONPATH=$PYTHONPATH:~/py
# cp -f ./System/Makefile.compiler.gfortran ./Makefile.compiler
  - modify the compiler flags to your liking (e.g. optimization, debug flags etc.)
- a serial run is configured by default, mpimod.f is linked to mpimod_ser.f.  If you want to run in MPI mode:
  # ln -sf ./mpimod_mpi.f ./mpimod.f
  - make sure to update the ./Makefile.compiler to something like: FC=mpif90
# make -j
- report any errors to: daan at flash.uchicago.edu


SETUP TIPS
==========
In order to create multiple builds (e.g. serial and MPI) from a single source, use lndir in the following way:

# mkdir -p supernu/src
# mkdir -p supernu/bin
# mkdir -p supernu/bin-mpi
# hg clone https://drrossum@bitbucket.org/drrossum/supernu supernu/src
# cd supernu/bin
# lndir ../src
# cd supernu/bin-mpi
# lndir ../src


REFERENCES
==========
Wollaeger, van Rossum et al. 2013
Densmoore 201x
Kasen et al. 200x
