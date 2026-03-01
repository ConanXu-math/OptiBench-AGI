import numpy as np
from scipy.optimize import minimize

def generate_synthetic_data(seed: int = 42) -> dict:
    np.random.seed(seed)
    d = 5  # dimension of x
    alpha = 0.1  # scalar parameter
    x_init = np.random.rand(d)  # initial guess
    
    def f(x):
        return np.sum(np.sin(x))  # example objective
    
    return {
        'x': x_init,
        'alpha': alpha,
        'f': f,
        'd': d
    }


def solve(data: dict) -> dict:
    f = data['f']
    x0 = data['x']
    result = minimize(f, x0, method='Nelder-Mead')
    return {
        'optimal_value': result.fun,
        'solution': result.x
    }

if __name__ == '__main__':
    data = generate_synthetic_data()
    print('Generated data:', data)
    solution = solve(data)
    print('Optimal value:', solution['optimal_value'])
    print('Solution:', solution['solution'])