import numpy as np
from numba import njit

n_samples = 1000000
n_bootstrap_samples = 1000

mu_a = 0.1

print('#mu_a\tsigma_e\t<|epistasis|>/<|fitness_effect|>')
for sigma_e in [0.01, 0.1, 1.]:
    @njit
    def epistasis_fitness_effect():
        f00 = np.random.normal(0., sigma_e)
        f01 = mu_a + np.random.normal(0., sigma_e)
        f10 = mu_a + np.random.normal(0., sigma_e)
        f11 = 2. * mu_a + np.random.normal(0., sigma_e)

        epistasis = np.abs((f10 - f00) - (f11 - f01))
        fitness_effect = np.abs(f10 - f00)

        return (epistasis, fitness_effect)

    # Calculate epistasis and fitness effects
    epistasis = np.zeros(n_samples)
    fitness_effect = np.zeros(n_samples)
    for sample in range(n_samples):
        r = epistasis_fitness_effect()
        epistasis[sample]      = r[0]
        fitness_effect[sample] = r[1]

    # Calculate bootstrap samples
    bootstrap_samples = np.zeros(n_bootstrap_samples)
    for sample in range(n_bootstrap_samples):
        indices = np.random.choice(n_samples, size=n_samples)
        bootstrap_samples[sample] = np.mean(epistasis[indices])/np.mean(fitness_effect[indices])

    print(mu_a, sigma_e, np.mean(bootstrap_samples), np.std(bootstrap_samples), sep='\t')
