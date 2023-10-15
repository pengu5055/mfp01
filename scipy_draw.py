import matplotlib.pyplot as plt
import numpy as np
from scipy import special

# scipy.special.airy(z) returns ndarrays ai, aip, bi, bip

# print(f"{a:.30f}")

x = np.linspace(-15, 2, 201)
ai, aip, bi, bip = special.airy(x)
plt.plot(x, ai, c="#83de5f", label="Ai(x)")
plt.plot(x, bi, c="#ec8cf5", label="Bi(x)", ls="--")


plt.title("Airyjevi funkciji")
plt.axhline(alpha=1, ls=":", c="#adadad")
plt.ylim(-0.6, 0.6)
plt.legend()
plt.show()
