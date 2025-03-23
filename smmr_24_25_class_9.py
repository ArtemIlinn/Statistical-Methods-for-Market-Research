import numpy as np
import matplotlib.pyplot as plt

# Markov chain simulation
S = [1, 2, 3]  # State space
P = np.array([[0.5, 0.25, 0.25],
             [0.3333, 0, 0.6667],
             [0.5, 0.5, 0]])  # Transition matrix

# Initialize chain with first state
X = [np.random.choice(S)]

# Generate remaining states
for t in range(1, 1000):
    prev_state = X[t-1]
    next_probs = P[prev_state-1]  # Adjust for 0-based index
    X.append(np.random.choice(S, p=next_probs))

# Plot first 50 states
plt.figure(figsize=(10, 4))
plt.plot(X[:50], drawstyle='steps-post', color='blue')
plt.xlabel('Time Step')
plt.ylabel('State')
plt.title('First 50 Steps of Markov Chain Simulation')
plt.yticks([1, 2, 3])
plt.grid(True)
plt.show()

# Transition count matrices initialization
P1_est = np.zeros((3, 3), dtype=int)
P2_est = np.zeros((3, 3), dtype=int)

# Count transitions for first half (t=2-500)
for t in range(1, 500):  # Python is 0-based
    prev = X[t-1] - 1  # Convert to 0-based index
    curr = X[t] - 1
    P1_est[prev, curr] += 1

# Count transitions for second half (t=501-1000)
for t in range(500, 1000):  # Python is 0-based
    prev = X[t-1] - 1  # Convert to 0-based index
    curr = X[t] - 1
    P2_est[prev, curr] += 1

print("First Half Transition Count Matrix:")
print(P1_est)
print("\nSecond Half Transition Count Matrix:")
print(P2_est)


# Convert counts to probabilities
P1_prob = P1_est / P1_est.sum(axis=1, keepdims=True)
P2_prob = P2_est / P2_est.sum(axis=1, keepdims=True)

print("\nFirst Half Transition Probabilities:")
print(np.round(P1_prob, 4))
print("\nSecond Half Transition Probabilities:")
print(np.round(P2_prob, 4))