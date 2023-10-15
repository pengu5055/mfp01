import matplotlib.pyplot as plt
import numpy as np
from scipy import special

# scipy.special.airy(z) returns ndarrays ai, aip, bi, bip

# print(f"{a:.30f}")


def f(array, n):  # Make recursive pls
    result = []
    for x in np.nditer(array):
        c = 0
        for k in range(0, n+1):
            step = special.gamma(1/3 + k)/special.gamma(1/3) * (3**k * x**(3*k))/(special.factorial(3*k))
            c += step
        result.append(c)
    return np.array(result)


def g(array, n):
    result = []
    for x in np.nditer(array):
        c = 0
        for k in range(0, n+1):
            step = special.gamma(2/3 + k)/special.gamma(2/3) * (3**k * x**(3*k + 1))/(special.factorial(3*k + 1))
            c += step
        result.append(c)
    return np.array(result)


def Ai(x, alpha, beta, n):
    return alpha*f(x, n) - beta*g(x, n)


def Bi(x, alpha, beta, n):
    return np.sqrt(3)*(alpha*f(x, n) + beta*g(x, n))


a = 0.355028053887817239
b = 0.258819403792806798
x_lin = np.linspace(-15, 5, 201)
ai, aip, bi, bip = special.airy(x_lin)
plt.plot(x_lin, ai, c="#83de5f", label="Ai(x)")
plt.plot(x_lin, bi, c="#ec8cf5", label="Bi(x)", ls="--")
plt.plot(x_lin, Ai(x_lin, a, b, 40), c="#ba2929")
plt.plot(x_lin, Bi(x_lin, a, b, 40), c="#60d4f7")

plt.title("Airyjevi funkciji")
plt.axhline(alpha=1, ls=":", c="#adadad")
plt.ylim(-0.6, 0.6)
plt.legend()
plt.show()
