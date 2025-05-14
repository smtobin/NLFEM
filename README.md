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
## Final Project
Working implementation of finite strain plasticity with isotropic hardening in 2D. 

Parameters: 1x1 2D square element, E=12000, nu=0.3, yield stress=100, H_bar_prime=1000, 10 load steps.

Two cases were tested:

### Case 1 - Uniaxial Stress
To run (from build directory):

```
./NLFEM ../input/final_case1.txt
```

Left nodes have x-displacement fixed, prescribed x-displacement on right nodes of 0.04. Bottom nodes have y-displacement fixed, y-displacement of top nodes free.

![image](https://github.com/user-attachments/assets/688e11c4-63c1-43e2-b2f7-aa111c4192d6)

Stresses and alpha for integration point 3:

![Screenshot from 2025-05-14 10-38-07](https://github.com/user-attachments/assets/8dc0c1a4-33e6-4241-868f-8b3c80dd8d7c)

### Case 2 - Right-side loading
To run (from build directory):

```
./NLFEM ../input/final_case2.txt
```

Left nodes are completely fixed, prescribed x-displacement on right nodes of 0.04. Y-displacement of right nodes free.
![image](https://github.com/user-attachments/assets/5c5fbc8f-f5fe-4ffe-9a99-8349c2ac9003)

Stresses and alpha for integration point 3:

![image](https://github.com/user-attachments/assets/43686d89-a3dc-4894-a78b-3fb37a256ab5)

### Why I think my results are correct

In the first case, the results agree with those posted on Canvas, and my code achieves these results with quadratic convergence.

In the second case, the results agree with my intuition (y-displacements are symmetric, the top-right node moves down and the bottom-right node moves up).

In both cases, I converge to the same results regardless of the number of load steps used.


## HW6
The repository for HW6: https://github.com/smtobin/NLFEM/tree/5f17ae68db1e84209f567faeb99ace00bae1dd4b

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
