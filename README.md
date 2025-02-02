# Code for ME 661 - Nonlinear FEM
## Linear FEM Assignment
To build the code (assuming Linux):
```
git clone git@github.com:smtobin/NLFEM.git
cd NLFEM
mkdir build
cd build
cmake ..
make
```
This will create an executable `LinearFEM` in the `build/` directory, which takes no inputs (everything is hard-coded for now).
```
./LinearFEM
```
This should print out all the correct info at the integration points for the Linear FEM assignment.
