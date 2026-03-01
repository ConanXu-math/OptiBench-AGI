import numpy as np
from scipy.optimize import minimize
import warnings

# Generate synthetic test instance
def generate_synthetic_data(seed: int = 42) -> dict:
    np.random.seed(seed)
    n, m = 5, 3  # dimensions for z and constraints

    # Define action set D (for x variable)
    D = np.random.randn(10, n)  # |D| = 10 actions
    bar_x_D = np.linalg.norm(D, axis=1).max()  # maximum Euclidean norm of points in D

    # Define feasibility region C
    # Use quadratic constraints for simplicity, ensuring Slater point exists
    Q = np.eye(n)  # identity matrix
    A = np.random.randn(m, n)  # constraint matrix
    b = np.ones(m)  # right-hand side
    bar_z = np.zeros(n)  # Slater point: satisfies A@bar_z < b
    
    # Ensure bar_z is strictly feasible
    while np.any(A @ bar_z >= b):
        bar_z += 0.1 * np.random.randn(n)
    
    # Define curvature constants μ_f, μ_g
    mu_f = 1.0  # curvature constant for f
    mu_g = np.ones(m)  # vector of curvature constants
    
    # Define upper bound on Lagrange multipliers
    bar_lambda = 10.0
    
    # Parameters α (dual step size), β (primal step size)
    alpha = 0.1
    beta = 0.5

    return {
        'z': np.zeros(n),  # initial z
        'x': D,            # action set
        'lambda': np.zeros(m),
        'beta': beta,
        'alpha': alpha,
        'mu_f': mu_f,
        'mu_g': mu_g,
        'bar_z': bar_z,
        'bar_lambda': bar_lambda,
        'bar_x_D': bar_x_D,
        'A': A,             # constraint matrix
        'b': b,             # right-hand side
        'm': m,             # number of constraints
        'n': n              # dimension of z
    }

# Solve the optimization problem
def solve(data: dict) -> dict:
    
    z_initial = data['z']
    A = data['A']
    b = data['b']
    n = data['n']
    m = data['m']
    alpha = data['alpha']
    beta = data['beta']
    mu_f = data['mu_f']
    mu_g = data['mu_g']
    bar_lambda = data['bar_lambda']
    bar_z = data['bar_z']
    
    def objective(z):
        # Simple quadratic function f(z) = ||z||²/2 + linear term to make it non-trivial
        q = 0.5 * np.dot(z, z)  # convex part
        linear_term = -np.sum(z[:n])  # additive linear component
        return q + linear_term

    def constraint_function(x):
        # Vector of constraint values: g(z) <= 0
        return A @ x  # returns an array of shape (m, ) → needs to be <= 0

    # Define constrained minimization problem
    res = minimize(
        fun=objective,
        x0=z_initial,
        method='SLSQP',
        constraints={'type': 'ineq', 'fun': lambda x: -constraint_function(x)},
        options={'disp': False, 'ftol': 1e-6, 'eps': 1e-9}
    )
    
    if not res.success:
        warnings.warn(f"Solver failed: {res.message}")
        # Fallback: use feasible point as solution
        solution = bar_z
        optimal_value = objective(solution)
    else:
        solution = res.x
        optimal_value = res.fun

    # Prepare output with minimal structure required
    result = {
        'optimal_value': float(optimal_value),
        'solution': solution.tolist(),
        'success': res.success,
        'message': res.message
    }
    
    return result

if __name__ == '__main__':
    # Generate test data
    data = generate_synthetic_data()
    print("Generated Data:")
    for key, val in data.items():
        print(f"{key}: shape={val.shape if isinstance(val, np.ndarray) else type(val)}")

    # Solve
    solution_dict = solve(data)
    
    print('\nSolution:')
    print(f"Optimal Value: {solution_dict['optimal_value']:.6f}")
    print(f"Solution Vector (z): {solution_dict['solution']}\n")