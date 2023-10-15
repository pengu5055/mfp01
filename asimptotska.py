import matplotlib.pyplot as plt
import numpy as np
from scipy import special

epsilon = 10**-10

def u_s(s):  # Tista zadeva k jo L, P, Q potem uporabljajo
    return special.gamma(3* s + 1/2)/(special.factorial(s) * 54**s * special.gamma(s + 1/2))

def xi(x):
    return 2/3 * np.abs(x)**(3/2)


def L(array, n):
    result = []
    for x in np.nditer(array):
        c = 0  # Counter
        for s in range(0, n+1):
            step = u_s(s)/x**s
            c += step
        result.append(c)
    return np.array(result)


def P(array, n):
    result = []
    for x in np.nditer(array):
        c = 0  # Counter
        for s in range(0, n + 1):
            step = (-1)**s * u_s(2*s)/x**(2*s)
            c += step
        result.append(c)
    return np.array(result)


def Q(array, n):
    result = []
    for x in np.nditer(array):
        c = 0  # Counter
        for s in range(0, n + 1):
            step = (-1)**s * u_s(2*s + 1) / x ** (2*s + 1)
            c += step
        result.append(c)
    return np.array(result)


def Ai_pos(array, n):
    return np.exp(-xi(array))/(2*np.sqrt(np.pi) * array**(1/4)) * L(-xi(array), n)


def Bi_pos(array, n):
    return np.exp(xi(array)) / (np.sqrt(np.pi) * array ** (1 / 4)) * L(xi(array), n)


def Ai_neg(array, n):
    return 1/(np.sqrt(np.pi)*(-array)**(1/4)) * (np.sin(xi(array) - np.pi/4)*Q(xi(array), n) + np.cos(xi(array) - np.pi/4)*P(xi(array), n))


def Bi_neg(array, n):
    return 1/(np.sqrt(np.pi)*(-array)**(1/4)) * (-np.sin(xi(array) - np.pi/4)*P(xi(array), n) + np.cos(xi(array) - np.pi/4)*Q(xi(array), n))

a = 0.355028053887817239
b = 0.258819403792806798
x_lin = np.linspace(1, 10, 201)
x_lin2 = np.linspace(-200, -2, 20001)
ai, aip, bi, bip = special.airy(x_lin)
ai2, aip2, bi2, bip2 = special.airy(x_lin2)
plt.plot(x_lin, ai, c="#83de5f", label="Ai(x)")
plt.plot(x_lin, bi, c="#ec8cf5", label="Bi(x)", ls="--")
plt.plot(x_lin2, ai2, c="#83de5f")
plt.plot(x_lin2, bi2, c="#ec8cf5", ls="--")
plt.plot(x_lin, Ai_pos(x_lin, 10), c="#ba2929")
plt.plot(x_lin, Bi_pos(x_lin, 15), c="#60d4f7")
plt.plot(x_lin2, Ai_neg(x_lin2, 10), c="#ba2929")
plt.plot(x_lin2, Bi_neg(x_lin2, 10), c="#60d4f7")

plt.title("Airyjevi funkciji")
plt.axhline(alpha=1, ls=":", c="#adadad")
plt.ylim(-0.6, 0.6)
plt.legend()
plt.show()