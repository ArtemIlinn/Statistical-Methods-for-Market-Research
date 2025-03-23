import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from scipy.optimize import minimize
from scipy.special import expit, logit

# Disable scientific notation
pd.set_option('display.float_format', lambda x: '%.5f' % x)

# Problem 1: Simple Logistic Regression
# =====================================

# Data simulation
np.random.seed(42)
n_obs = 1000
beta_0 = 2
beta_1 = -3

epsilon = stats.logistic.rvs(size=n_obs)
x = stats.bernoulli.rvs(0.75, size=n_obs)
y = (beta_0 + beta_1*x + epsilon > 0).astype(int)
df = pd.DataFrame({'y': y, 'x': x})

# Logit model estimation
X = sm.add_constant(df[['x']])
logit_model = sm.Logit(df['y'], X).fit()
print("Problem 1 Coefficients:")
print(logit_model.params)

# Analytical solution
n00 = ((y == 0) & (x == 0)).sum()
n01 = ((y == 0) & (x == 1)).sum()
n10 = ((y == 1) & (x == 0)).sum()
n11 = ((y == 1) & (x == 1)).sum()

analytical_beta0 = np.log(n10/n00)
analytical_beta1 = np.log((n11*n00)/(n01*n10))

print("\nAnalytical Solutions:")
print(f"Beta0: {analytical_beta0:.4f}")
print(f"Beta1: {analytical_beta1:.4f}")

# Problem 2: Logistic Regression with Marginal Effects
# ====================================================

# Data simulation
np.random.seed(42)
n_obs = 1000
beta_0 = 2
beta_1 = -3
beta_2 = 3

epsilon = stats.logistic.rvs(size=n_obs)
x1 = stats.expon.rvs(size=n_obs)
x2 = stats.bernoulli.rvs(0.75, size=n_obs)
y = (beta_0 + beta_1*x1 + beta_2*x2 + epsilon > 0).astype(int)
df = pd.DataFrame({'y': y, 'x1': x1, 'x2': x2})

# Logit model estimation
X = sm.add_constant(df[['x1', 'x2']])
logit_model2 = sm.Logit(df['y'], X).fit()
print("\nProblem 2 Coefficients:")
print(logit_model2.params)

# Marginal effects (average marginal effects)
predicted_probs = logit_model2.predict(X)
marginal_effects = pd.DataFrame({
    'x1': logit_model2.params['x1'] * predicted_probs * (1 - predicted_probs),
    'x2': logit_model2.params['x2'] * predicted_probs * (1 - predicted_probs)
})

print("\nMarginal Effects (Average):")
print(marginal_effects.mean())

# Custom Log-Likelihood Function
# =============================

def log_likelihood(beta, y, x1, x2):
    Xb = beta[0] + beta[1]*x1 + beta[2]*x2
    log_lik = np.sum(y * np.log(expit(Xb))) + np.sum((1-y) * np.log(1 - expit(Xb)))
    return -log_lik  # Negative for minimization

# Numerical optimization
result = minimize(log_likelihood,
                  x0=np.zeros(3),
                  args=(df['y'], df['x1'], df['x2']),
                  method='Nelder-Mead')

print("\nOptimization Results:")
print("Coefficients:", result.x)
print("Log-Likelihood:", -result.fun)

# Calculate Hessian-based standard errors
hessian = result.hess_inv(np.eye(3)) if result.hess_inv is not None else np.zeros((3,3))
std_errors = np.sqrt(np.diag(hessian))
print("\nStandard Errors:", std_errors)