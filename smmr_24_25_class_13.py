import numpy as np
from scipy.stats import kstest, pareto, genextreme, norm, gumbel_r
import matplotlib.pyplot as plt

# Problem 1
# -------------------------------------------------
n_obs = 1000
n_sim = 1000
theta = 3

# Simulate Pareto maxima
pareto_sam_max = np.empty(n_sim)
for m in range(n_sim):
    pareto_sam_aux = theta * (pareto.rvs(1, size=n_obs) + 1)  # Pareto(II) implementation
    pareto_sam_max[m] = np.max(pareto_sam_aux)

c_n = theta * n_obs
alpha = 1

# Generate Fréchet samples (genextreme with c=-1/alpha corresponds to Fréchet)
frechet_sam = genextreme.rvs(c=-1/alpha, loc=0, scale=1, size=n_obs)

# KS test
ks_stat, p_value = kstest(pareto_sam_max / c_n, frechet_sam)
print(f"Problem 1 KS Test: statistic={ks_stat:.4f}, p-value={p_value:.4f}")

# Plot histograms
plt.hist(pareto_sam_max / c_n, bins=100, density=True, alpha=0.5, label='Normalized Pareto Maxima')
plt.hist(frechet_sam, bins=100, density=True, alpha=0.5, label='Fréchet Samples')
plt.legend()
plt.title('Problem 1 Distribution Comparison')
plt.show()

# Problem 2
# -------------------------------------------------
n_obs = 100000
n_sim = 1000

# Simulate Normal maxima
norm_sam_max = np.empty(n_sim)
for m in range(n_sim):
    norm_sam_aux = norm.rvs(size=n_obs)
    norm_sam_max[m] = np.max(norm_sam_aux)

# Calculate normalization constants
log_term = np.log(n_obs)
c_n = 1 / np.sqrt(2 * log_term)
d_n = np.sqrt(2 * log_term) - (np.log(4 * np.pi) + np.log(log_term)) / (2 * np.sqrt(2 * log_term))

# Generate Gumbel samples
gumbel_sam = gumbel_r.rvs(size=n_obs)

# KS test
ks_stat, p_value = kstest((norm_sam_max - d_n) / c_n, gumbel_sam)
print(f"Problem 2 KS Test: statistic={ks_stat:.4f}, p-value={p_value:.4f}")

# Plot histograms
plt.hist((norm_sam_max - d_n) / c_n, bins=100, density=True, alpha=0.5, label='Normalized Normal Maxima')
plt.hist(gumbel_sam, bins=100, density=True, alpha=0.5, label='Gumbel Samples')
plt.hist(norm_sam_max, bins=10, density=True, alpha=0.3, label='Original Maxima')
plt.legend()
plt.title('Problem 2 Distribution Comparison')
plt.show()