import numpy as np
from scipy.stats import norm
from scipy.optimize import minimize
import pygad

# ======================
# TASK 1: Simple MSM Model
# ======================

# Simulation parameters
n_obs = 2000
np.random.seed(42)

# Markov chain simulation
S = [0, 1]  # States 0 and 1 (Python uses 0-based indexing)
P = np.array([[0.95, 0.05],
              [0.15, 0.85]])

M = np.zeros(n_obs, dtype=int)
M[0] = np.random.choice(S)
for i in range(1, n_obs):
    M[i] = np.random.choice(S, p=P[M[i-1]])

# Sample simulation
mu = [2, -2]
sigma = [1, 2]
X = np.zeros(n_obs)
for i in range(n_obs):
    state = M[i]
    X[i] = np.random.normal(mu[state], sigma[state])

# Log-likelihood function
def logL_SMM(par):
    mu1, mu2, sigma1, sigma2, p11, p22 = par
    P_t = np.array([[p11, 1-p11],
                    [1-p22, p22]])
    
    xi = np.array([0.5, 0.5])
    log_lik = 0
    
    for x in X:
        f = np.array([
            norm.pdf(x, mu1, sigma1),
            norm.pdf(x, mu2, sigma2)
        ])
        xi_next = P_t.T @ xi
        L_t = f @ xi_next
        log_lik += np.log(L_t)
        xi = (f * xi_next) / L_t
    
    return -log_lik  # Negative for minimization

# Constraints and optimization
bounds = [
    (None, None),   # mu1
    (None, None),   # mu2
    (0.01, None),   # sigma1
    (0.01, None),   # sigma2
    (0.01, 0.99),   # p11
    (0.01, 0.99)    # p22
]

initial_guess = [0, 0, 1, 1, 0.5, 0.5]
result = minimize(logL_SMM, initial_guess, 
                 method='SLSQP', bounds=bounds)

# Results
print("Estimated parameters:", result.x)
P_est = np.array([
    [result.x[4], 1-result.x[4]],
    [1-result.x[5], result.x[5]]
])

# Stationary distribution
pi = np.array([0.5, 0.5])
for _ in range(100):
    pi = pi @ P_est
print("Stationary distribution:", pi)

# ======================
# TASK 2: MSM for AR(1)
# ======================

# Simulation parameters
T_obs = 1000
np.random.seed(42)

# Markov chain simulation
M_ar = np.zeros(T_obs, dtype=int)
M_ar[0] = np.random.choice(S)
for t in range(1, T_obs):
    M_ar[t] = np.random.choice(S, p=P[M_ar[t-1]])

# AR(1) simulation
alpha = [0.1, 0.9]
beta = [0.8, 0.3]
sigma_ar = [2, 5]
X_ar = np.zeros(T_obs)
for t in range(1, T_obs):
    state = M_ar[t]
    X_ar[t] = alpha[state] + beta[state]*X_ar[t-1] + np.random.normal(0, sigma_ar[state])

# Log-likelihood function
def fitness_func(ga_instance, solution, solution_idx):
    alpha1, alpha2, beta1, beta2, sigma1, sigma2, p11, p22 = solution
    P_t = np.array([[p11, 1-p11],
                    [1-p22, p22]])
    
    xi = np.array([0.5, 0.5])
    log_lik = 0
    
    for t in range(T_obs-1):
        x_next = X_ar[t+1]
        x_prev = X_ar[t]
        
        f = np.array([
            norm.pdf(x_next, alpha1 + beta1*x_prev, sigma1),
            norm.pdf(x_next, alpha2 + beta2*x_prev, sigma2)
        ])
        
        xi_next = P_t.T @ xi
        L_t = f @ xi_next
        log_lik += np.log(L_t)
        xi = (f * xi_next) / L_t
    
    return log_lik

# Genetic algorithm setup
ga_instance = pygad.GA(
    num_generations=200,
    num_parents_mating=10,
    fitness_func=fitness_func,
    sol_per_pop=150,
    num_genes=8,
    gene_space=[
        (-10, 10), (-10, 10),   # alphas
        (-0.99, 0.99), (-0.99, 0.99),  # betas
        (0.01, 10), (0.01, 10),  # sigmas
        (0.01, 0.99), (0.01, 0.99)  # p11, p22
    ],
    mutation_probability=0.1,
    parent_selection_type="sss",
    crossover_type="single_point",
    mutation_type="random",
    keep_parents=2,
    stop_criteria="saturate_100"
)

# Run optimization
ga_instance.run()

# Results
solution = ga_instance.best_solution()
print("Best parameters:", solution[0])
print("Log-likelihood:", solution[1])