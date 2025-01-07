import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt

def elliptic_solver(Mx, My):
    # Define the spatial domain bounds
    ax, bx = 0, 1
    ay, by = 0, 1

    # Compute grid spacing
    hx = (bx - ax) / Mx
    hy = (by - ay) / My

    # Generate grid points
    x = np.linspace(ax, bx, Mx + 1)
    y = np.linspace(ay, by, My + 1)
    xx, yy = np.meshgrid(x, y)

    # Linear index function
    def myind(i, j):
        return (j - 1) * (Mx + 1) + i

    # Initialize the sparse matrix, RHS vector, and solution vector
    Nall = (Mx + 1) * (My + 1)
    A = lil_matrix((Nall, Nall))
    F = np.zeros(Nall)
    U = np.zeros(Nall)
    Uexact = np.zeros(Nall)

    # Define the right-hand side function
    def rhsfun(x, y):
        return 2  # Example forcing function

    # Define the exact solution for validation
    def exactfun(x, y):
        return x * (1 - x) / 2 + y * (1 - y) / 2

    # Build the sparse matrix and RHS vector
    interior_nodes = []
    for i in range(1, Mx):
        for j in range(1, My):
            me = myind(i, j)
            interior_nodes.append(me)
            F[me] = rhsfun(x[i], y[j])
            Uexact[me] = exactfun(x[i], y[j])

            # Neighboring indices
            mel = myind(i - 1, j)  # left
            mer = myind(i + 1, j)  # right
            meb = myind(i, j - 1)  # bottom
            met = myind(i, j + 1)  # top

            A[me, me] = 2 / hx**2 + 2 / hy**2
            A[me, mel] = -1 / hx**2
            A[me, mer] = -1 / hx**2
            A[me, meb] = -1 / hy**2
            A[me, met] = -1 / hy**2

    # Apply boundary conditions
    for j in range(My + 1):  # Left and right boundaries
        A[myind(0, j), myind(0, j)] = 1  # Left boundary
        A[myind(Mx, j), myind(Mx, j)] = 1  # Right boundary

    for i in range(Mx + 1):  # Bottom and top boundaries
        A[myind(i, 0), myind(i, 0)] = 1  # Bottom boundary
        A[myind(i, My), myind(i, My)] = 1  # Top boundary

    # Solve the sparse linear system
    A = A.tocsr()  # Convert to compressed sparse row format
    U = spsolve(A, F)

    # Reshape the solution into 2D arrays for plotting
    Uplot = np.zeros((Mx + 1, My + 1))
    Uexactplot = np.zeros((Mx + 1, My + 1))
    for i in range(Mx + 1):
        for j in range(My + 1):
            me = myind(i, j)
            Uplot[i, j] = U[me]
            Uexactplot[i, j] = Uexact[me]

    # Plot the numerical solution
    plt.figure()
    ax1 = plt.axes(projection='3d')
    ax1.plot_surface(xx, yy, Uplot.T, cmap='viridis')
    plt.title('Numerical Solution')
    plt.show()

    # Calculate and print the error
    error = np.linalg.norm(Uexact - U, np.inf)
    print(f"Error = {error}")

# Example usage
elliptic_solver(10, 10)
