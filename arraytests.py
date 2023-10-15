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

x_lin = np.linspace(-15, 2, 201)
print(Ai(x_lin, a, b,10))
plt.plot(x_lin, Ai(x_lin, a, b, 10))
plt.show()
