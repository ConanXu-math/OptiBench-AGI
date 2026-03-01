import numpy as np
from scipy.optimize import minimize


def generate_synthetic_data(seed: int = 42) -> dict:
    np.random.seed(seed)
    n = 5  # dimension of x
    # Create a compact feasible set: a unit cube [0, 1]^n
    x_low = np.zeros(n)
    x_high = np.ones(n)
    x = np.random.uniform(x_low, x_high, size=n)  # feasible point inside K
    
    # Construct objective: f(x) = sum_{i=1}^{n} x_i^2 + sum_{i<j} x_i * x_j
    # This is a quadratic (degree 2) polynomial, convex if symmetric and positive definite
    # But note: the problem says polynomial of degree ≥ 2 — this is acceptable.
    # We construct it to be strictly convex.
    A = np.eye(n) + np.random.rand(n, n)  # ensure positive definite (convex quadratic)
    b = np.zeros(n)
    c = 0.0
    
    return {
        'n': n,
        'x': x,
        'A': A,  # Hessian matrix for quadratic part
        'b': b,
        'c': c,
        'K_bounds': (x_low, x_high),
        'K_dim': n
    }


def solve(data: dict) -> dict:
    n = data['n']
    x0 = data['x']
    A = data['A']
    b = data['b']
    c = data['c']
    x_low, x_high = data['K_bounds']
    
    def objective(x):
        # quadratic objective: (1/2) x^T A x + b^T x + c
        # We use 1/2 to make derivative nice for minimization
        return 0.5 * np.dot(x.T, np.dot(A, x)) + np.dot(b, x) + c
    
    def constraint_fn(x):
        # Enforce bounds: x >= x_low, x <= x_high
        return np.concatenate([x - x_low, x_high - x])
    
    constraints = ({
        'type': 'ineq',
        'fun': lambda x: constraint_fn(x)
    })
    
    result = minimize(
        fun=objective,
        x0=x0,
        method='SLSQP',
        bounds=list(zip(x_low, x_high)),
        constraints=constraints,
        options={'disp': True}
    )
    
    return {
        'optimal_value': result.fun,
        'solution': result.x,
        'success': result.success,
        'message': result.message
    }

if __name__ == '__main__':
    data = generate_synthetic_data()
    solution = solve(data)
    print('Optimal Value:', solution['optimal_value'])
    print('Solution Vector:', solution['solution'])
    print('Success:', solution['success'])