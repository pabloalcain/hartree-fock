#Hartree-Fock solution of atoms

This software is used in Quantum Chemistry course in University of
Buenos Aires. For more reference, go to the
[webpage](http://materias.df.uba.ar/e3a2016c1/) [in Spanish!]. It
ports a Hartree-Fock solution for central potentials made by 
[Walter Johnson][1] for the Notre Dame University

[1]: http://www3.nd.edu/~johnson/

DISCLAIMER: The files here do not follow any guideline respecting to
the usual python module distribution. It is done so in purpose, so
students *know* what the sources are.


## Installation
After downloading the package, go to `src/` folder and run `make
library`

```
cd src/
make library
```

## Running

In the root folder of the project, there is an example called
`krypton.py` that calculates Hartree-Fock solution for a neutral Kr
atom. You can (*and should*) take a look at its source code and
modify it properly.

## Help
The file `hartree-fock.py` has the help embedded in the docstrings

## Another way
In order to keep the previous format working, you can also build the old
(and less user friendly) executable file with:

```
cd src/
make executable
```

After that, you can run the file `./johnson.e` and work with the
executable file directly. Keep in mind that if you do that, you'll
have to manually extract all info like the wavefunction and the
energy. We insist on *not using* this file, except it is strictly
necessary.

## Note

The `johnson.for` file had no license information and no reference to
link to, so I left the header and comments. If this violates any
copyright, let me know and I will take this down. Please, if you
share, keep the copyright information on the original (and the hardest
part of the code!) `johnson.for` file.
