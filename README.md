# Code for ME 661 - Nonlinear FEM
## Building and Running the Code
To build the code (assuming Linux):
```
git clone git@github.com:smtobin/NLFEM.git
cd NLFEM
mkdir build
cd build
cmake ..
make
```
This will create an executable `NLFEM` in the `build/` directory, which takes the path to an input file as input.
```
./NLFEM ../input/linear_fem.txt
```
## Linear FEM Assignment
The Linear FEM assignment can be ran with:
```
./NLFEM ../input/linear_fem.txt
```

which will print out the correct info at the integration points of the single element.

## HW3
HW3 assignment consists of 3 parts (a), (b) and (c) which each have their own input files:
```
./NLFEM ../input/hw3a.txt
./NLFEM ../input/hw3b.txt
./NLFEM ../input/hw3c.txt
```
Each input file should be set up to satisfy the question requirements, and prints out info at each integration point of the single element.

## Input File Format
The input files specify the material, mesh nodes, mesh elements, prescribed displacements and prescribed forces. An example input file is below:
```
# material - choose from "PlaneStrain" or "HW3"
m PlaneStrain

# nodes - index, x, y
n 0 0.5 0.5
n 1 -0.5 0.5
n 2 -0.5 -0.5
n 3 0.5 -0.5

# elements - n1, n2, n3, n4
e 0 1 2 3

# displacement BC - node index, axis (x or y), displacement
d 0 x 0.01
d 1 x 0
d 2 x 0
d 3 x 0.01
d 3 y 0

# force BC - node index, axis (x or y), force
# example: f 0 x 10
# f 0 x 1.655
# f 1 x -1.655
# f 2 x -1.655
# f 3 x 1.655
```
