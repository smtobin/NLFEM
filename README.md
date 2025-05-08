# Code for ME 661 - Nonlinear FEM
## Building and Running the Code
First, install GNUPlot with
```
sudo apt-get install gnuplot
```

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
./NLFEM ../input/hw6.txt
```
## HW6
To run the radial return algorithm for Q1, run
```
./RadialReturn
```
from the build folder. The generated graphs are shown below:

![image](https://github.com/user-attachments/assets/79617f7e-621a-4245-afc4-572ab809af8b)


To run the full FEM solver for Q2, run
```
./NLFEM ../input/hw6.txt
```
from the build folder. The generated graphs are shown below:

![image](https://github.com/user-attachments/assets/191ac4f0-1606-4107-b4ee-1f4de8478db7)


The applied boundary conditions are a pin joint at the bottom left node, setting the x displacement on the upper left node, and prescribing 0.02 x displacement on the right two nodes. The number of load steps is hard-coded to be 100, and the time step is hard-coded to be 0.01.

## Midterm Assignment
The repository for Midterm Assignment: https://github.com/smtobin/NLFEM/tree/28410df3ba68217287d24c6361fe2258d310452e

## Linear FEM Assignment
The repository for Linear FEM: https://github.com/smtobin/NLFEM/tree/3b38613a801e5e5c18c425d65bef1809654b406a

## HW3
The repository for HW3: https://github.com/smtobin/NLFEM/tree/d4e3f46f18a2723ba919e170b7eca36f22bb4d5b

## Input File Format
The input files specify the material, mesh nodes, mesh elements, prescribed displacements and prescribed forces. An example input file is below:
```
# material - only option is "Midterm"
m Midterm

# number of load steps
l 10

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
