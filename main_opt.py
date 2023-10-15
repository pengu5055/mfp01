import matplotlib.pyplot as plt
import numpy as np
from scipy import special
import mpmath

# Set percision
epsilon = 10**-10
mpmath.mp.dps = 10


def u_s(s):  # Tista zadeva k jo L, P, Q potem uporabljajo
    return mpmath.gamma(3 * s + 1/2)/(mpmath.fac(s) * 54**s * mpmath.gamma(s + 1/2))


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
    return mpmath.exp(-xi(array))/(2*mpmath.sqrt(np.pi) * array**(1/4)) * L(-xi(array), n)


def Bi_pos(array, n):
    return mpmath.exp(xi(array)) / (mpmath.sqrt(np.pi) * array ** (1 / 4)) * L(xi(array), n)


def Ai_neg(array, n):
    return 1/(mpmath.sqrt(mpmath.pi)*(-array)**(1/4)) * (mpmath.sin(xi(array) - mpmath.pi/4)*Q(xi(array), n) + mpmath.cos(xi(array) - mpmath.pi/4)*P(xi(array), n))


def Bi_neg(array, n):
    return 1/(mpmath.sqrt(mpmath.pi)*(-array)**(1/4)) * (-mpmath.sin(xi(array) - mpmath.pi/4)*P(xi(array), n) + mpmath.cos(xi(array) - mpmath.pi/4)*Q(xi(array), n))


def f(array, n):  # Make recursive pls
    result = []
    for x in np.nditer(array):
        c = 0
        for k in range(0, n+1):
            step = mpmath.gamma(1/3 + k)/mpmath.gamma(1/3) * (3**k * x**(3*k))/(mpmath.fac(3*k))
            c += step
        result.append(c)
    return np.array(result)


def g(array, n):
    result = []
    for x in np.nditer(array):
        c = 0
        for k in range(0, n+1):
            step = mpmath.gamma(2/3 + k)/mpmath.gamma(2/3) * (3**k * x**(3*k + 1))/(mpmath.fac(3*k + 1))
            c += step
        result.append(c)
    return np.array(result)


def Ai(x, alpha, beta, n):
    return alpha*f(x, n) - beta*g(x, n)


def Bi(x, alpha, beta, n):
    return mpmath.sqrt(3)*(alpha*f(x, n) + beta*g(x, n))


a = 0.355028053887817239
b = 0.258819403792806798
x_lin1 = np.linspace(-20, -5, 201)
x_lin2 = np.linspace(-5, 3, 201)
x_lin3 = np.linspace(3, 10, 201)

ai1, aip1, bi1, bip1 = special.airy(x_lin1)
ai2, aip2, bi2, bip2 = special.airy(x_lin2)
ai3, aip3, bi3, bip3 = special.airy(x_lin3)

plt.plot(x_lin1, ai1, c="#83de5f", label="Ai(x)")
plt.plot(x_lin1, bi1, c="#ec8cf5", label="Bi(x)", ls="--")
plt.plot(x_lin2, ai2, c="#83de5f")
plt.plot(x_lin2, bi2, c="#ec8cf5", ls="--")
plt.plot(x_lin3, ai3, c="#83de5f")
plt.plot(x_lin3, bi3, c="#ec8cf5", ls="--")

plt.plot(x_lin1, np.array(Ai_neg(x_lin1, 10)), c="#ba2929")
plt.plot(x_lin1, Bi_neg(x_lin1, 10), c="#60d4f7")

plt.plot(x_lin2, Ai(x_lin2, a, b, 25), c="#ba2929")
plt.plot(x_lin2, Bi(x_lin2, a, b, 25), c="#60d4f7")

plt.plot(x_lin3, Ai_pos(x_lin3, 10), c="#ba2929")
plt.plot(x_lin3, Bi_pos(x_lin3, 10), c="#60d4f7")

plt.title("Airyjevi funkciji")
plt.axhline(alpha=1, ls=":", c="#adadad")
plt.ylim(-0.6, 0.6)
plt.legend()
plt.show()