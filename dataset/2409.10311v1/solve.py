import numpy as np
from scipy.optimize import minimize


def generate_synthetic_data(seed: int = 42) -> dict:
    np.random.seed(seed)
    n = 10  # size of H
    m = 5   # size of G
    
    # Linear operator L (random matrix from R^{m x n})
    L = np.random.randn(m, n)
    
    # Define simple quadratic f and convex g for test instance
    def f(x):
        return np.sum(x**2) / 2
    
    def g(y):
        return np.sum(np.abs(y))  # l1 norm
    
    # Initialize variables with zeros
    x = np.zeros(n)
    y_hat = np.zeros(m)
    z_hat = np.zeros(m)
    
    # Set parameters
    gamma = 0.1
    tau = 0.8
    sigma = 1.0
    epsilon_k = 0.01
    
    # Return all data as dictionary
    return {
        'L': L,
        'f': f,
        'g': g,
        'gamma': gamma,
        'tau': tau,
        'sigma': sigma,
        'epsilon_k': epsilon_k,
        'x': x,
        'y_hat': y_hat,
        'z_hat': z_hat
    }


def solve(data: dict) -> dict:
    L = data['L']
    f = data['f']
    g = data['g']
    gamma = data['gamma']
    tau = data['tau']
    sigma = data['sigma']
    epsilon_k = data['epsilon_k']
    x_init = data['x']
    y_hat_init = data['y_hat']
    z_hat_init = data['z_hat']
    
    # Because this is an iterative algorithm, we need to do multiple steps
    x = x_init.copy()
    y_hat = y_hat_init.copy()
    z_hat = z_hat_init.copy()
    
    # We'll simulate one step since the problem is iterative and not a single optimization problem
    # In full ADMM-style implementation, you'd loop until convergence
    
    # Step 1: update x_k by solving:
    # minimize_x: f(x) + <z_hat | Lx - y_hat> + γ/2 ||Lx - y_hat||^2
    
    # This is equivalent to:
    # f(x) + (1/γ)[(Lx - y_hat) - z_hat] * [Lx - y_hat] + \text{const}
    
    # Simplify inner term w.r.t. x
    def obj_x(x_flat):
        x = x_flat.reshape(-1)
        fx = f(x)
        proj_term = np.dot(L.T, (L @ x - y_hat))  # gives gradient of Lx-y_hat terms
        quad_term = np.sum((L @ x - y_hat)**2) * gamma / 2
        linear_term = np.dot(z_hat, L @ x - y_hat)
        
        return fx + linear_term + quad_term
    
    # Use scipy.optimize.minimize
    res_x = minimize(obj_x, x.flatten(), method='BFGS')
    x_new = res_x.x.reshape(-1)
    
    # Step 2: Get v_k as subgradient of g at current ~y_k
    # Once you compute it
    # Assume that we take ~y_k = y_hat (simplification), so use subgradient of g
    # Since g is L1-norm, subgradient is sign function
    tilde_y = y_hat
    v = np.sign(tilde_y)  # subgradient of sum abs == l1 norm
    
    # Step 3: Compute residual e_k = v - z_hat + gamma*(tilde_y - L*x_new)
    e_k = v - z_hat + gamma * (tilde_y - L @ x_new)
    
    # Step 4: Check feasibility condition
    left_side = np.sum(e_k**2) + 2*gamma*epsilon_k
    right_side = sigma**2 * np.minimum(
        gamma**2 * np.sum((L @ x_new - tilde_y)**2),
        np.sum((v - z_hat)**2)
    )
    
    # Generate new iterates
    z_next = z_hat + tau * gamma * (L @ x_new - tilde_y)
    y_next = (1 - tau) * y_hat + tau / gamma * (z_hat + gamma * L @ x_new - v)
    
    # Package results
    result = {
        'optimal_value': f(x_new),  # crude estimate; real optimum requires more iterations
        'solution': x_new,
        'iteration_info': {
            'x': x_new,
            'y_hat': y_next,
            'z_hat': z_next,
            'v': v,
            'e_k': e_k,
            'converged': left_side <= right_side  # check if constraint satisfied within tolerance
        }
    }
    
    return result

if __name__ == '__main__':
    data = generate_synthetic_data(seed=42)
    sol = solve(data)
    print("Optimal Value:", sol['optimal_value'])
    print("Solution:", sol['solution'])
    print("Converged:", sol['iteration_info']['converged'])