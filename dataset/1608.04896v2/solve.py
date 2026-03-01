import numpy as np
from scipy.optimize import minimize

def robin_eigenvalue_approx(shape_params, c1, c2, alpha):
    """Approximates lowest Robin eigenvalue λ₁^(α)(Ω_ext) for given shape parameters."""
    a, b = shape_params
    
    # Ensure positive axes
    if a <= 0 or b <= 0:
        return float('-inf')
    
    # Compute perimeter and area of ellipse (approximations)
    perimeter = 4 * a * np.arcsin(np.sqrt((a**2 + b**2) / (2 * a * b))) if a == b else 2 * np.pi * np.sqrt((a**2 + b**2)/2)
    area = np.pi * a * b
    
    # Penalization if constraints violated
    perim_penalty = max(0, abs(perimeter - c1)) * 1e6
    area_penalty = max(0, abs(area - c2)) * 1e6
    
    # Simple heuristic for eigenvalue
    gamma = 1.0 + 0.1 * a / b
    eigenval_heuristic = gamma / (a ** 2 + b ** 2) if (a*b) > 0 else float('-inf')
    
    # Return objective (to be maximized => negate for minimizer)
    total_cost = -eigenval_heuristic + perim_penalty + area_penalty
    return total_cost

def generate_synthetic_data(seed=42) -> dict:
    """Generate synthetic data with feasible constraints."""
    np.random.seed(seed)
    
    # Feasible values
    c1 = 10.0  # perimeter constraint
    c2 = 6.0   # area constraint
    alpha = -0.5  # negative Robin parameter
    
    # Initialize with feasible starting point
    a_guess = 2.5
    b_guess = c2 / (np.pi * a_guess)  # ensure π*a*b ≈ c2
    
    return {
        'c1': c1,
        'c2': c2,
        'alpha': alpha,
        'a0': a_guess,
        'b0': b_guess
    }

def solve(data: dict) -> dict:
    """Solve optimization problem using scipy.optimize.minimize."""
    c1, c2, alpha = data['c1'], data['c2'], data['alpha']
    
    x0 = [data['a0'], data['b0']]  # ellipse parameters [a, b]
    bounds = [(0.1, 10), (0.1, 10)]  # avoid singular shapes
    
    def objective(x):
        return robin_eigenvalue_approx(x, c1, c2, alpha)
    
    result = minimize(
        fun=objective,
        x0=x0,
        bounds=bounds,
        method='SLSQP',
        options={'disp': True, 'ftol': 1e-8}
    )
    
    if not result.success:
        raise ValueError(f"Optimization failed: {result.message}")
    
    optimal_value = -result.fun
    solution = result.x
    
    return {
        'optimal_value': optimal_value,
        'solution': solution,
        'status': result.success,
        'message': result.message
    }

if __name__ == '__main__':
    data = generate_synthetic_data()
    print("Generated Data:")
    for k, v in data.items():
        print(f"{k}: {v}")
    
    solution = solve(data)
    print(f"Optimization Result:\nOptimal Value: {solution['optimal_value']:.4f}\nSolution (a, b): {solution['solution']}\nSuccess: {solution['status']}"')