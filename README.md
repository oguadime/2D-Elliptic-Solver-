# 2D-Elliptic-Solver-
A Python implementation of a 2D elliptic solver using the finite difference method. The solver is designed for numerical experimentation and supports validation with exact solutions.

# Description 
The script elliptic_solver solves the 2D elliptic equation: -(u_{xx}+u_{yy}) = f(x,y) in the domain (x,y) in [0,1]^2. The solver uses the finite difference method to discretize the domain and computes the numerical solution. It also supports visualization and error analysis against an exact solution.

# Key Features
- Custom Grid: Supports user-defined grid sizes along the x and y axes.
- Boundary Conditions: Implements Dirichlet boundary conditions.
- Exact Solution Comparison: Computes the infinity norm of the error if the exact solution is available.
- Visualization: Generates 3D surface plots of the numerical solution.
- Matrix Assembly: Constructs the sparse matrix for the finite difference discretization.

# Usage 
Syntax: 
```matlab 
elliptic_solver(Mx, My)
```

Inputs 
- Mx: Number of grid points along the x-axis (integer).
- My: Number of grid points along the y-axis (integer).

Outputs
- A 3D plot of the numerical solution.
- The infinity norm of the error compared to the exact solution is printed to the console.

# Examples 
Example with a 10x10 Grid:
```matlab
elliptic_solver(10, 10)
```

# License 
This project is licensed under the MIT License - see the LICENSE file for details.
```
Feel free to adjust any part of this README to better fit your specific needs or preferences.
```
