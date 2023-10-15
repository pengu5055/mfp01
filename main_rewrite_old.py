import matplotlib.pyplot as plt
from mpmath import *
import numpy as np

mp.dps = 100

def f(x, n):
    c = mpf(0)  # Counter
    x = mpf(x)
    for k in range(0, n+1):
        print("f: " + str(k))
        step = rf(1/3, k) * mpf("3")**mpf(k) * mpf(x)**mpf(3*k)/fac(3*k)
        old = c
        c += step
        if old == c:
            return c
    return c


def g(x, n):
    c = mpf(0)  # Counter
    x = mpf(x)
    for k in range(0, n+1):
        print("g: " + str(k))
        step = rf(2/3, k) * mpf("3")**mpf(k) * mpf(x)**mpf(3*k + 1)/fac(3*k + 1)
        old = c
        c += step
        if old == c:
            return c
    return c


def u_s(s):  # Deluje samo za mpf in ne array
    return gamma(mpf(3*s) + mpf("0.5"))/(54**s * fac(s)*gamma(s + mpf("0.5")))


def xi(x):
    return 2/3 * np.abs(x)**(3/2)


def L(x, n):
    c = mpf(0)
    x = mpf(x)
    oldstep = mpf(10e100)  # Samo nekaj zelo velikega, da ima na zacetku za primerjat z neko neskoncnostjo
    for s in range(0, n+1):
        step = u_s(s)/x**s
        c += step

    return c


def P(x, n):
    c = mpf(0)
    x = mpf(x)
    oldstep = mpf(10e100)
    for s in range(0, n+1):
        print("P: " + str(s))
        step = (-1)**s * u_s(2*s)/mpf(x)**(2*s)
        c += step
        #if abs(oldstep) <= abs(c):
        #    return c
        #else:
        #    oldstep = step
    return c


def Q(x, n):
    c = mpf(0)
    x = mpf(x)
    oldstep = mpf(10e100)
    for s in range(0, n+1):
        print("Q: " + str(s))
        step = (-1)**s * u_s(2*s + 1) / mpf(x) ** (2*s + 1)
        c += step
        #if abs(oldstep) <= abs(c):
        #    return c
        #else:
        #    oldstep = step
    return c


def arrayize(array, function, *args):
    return np.array([function(x, *args) for x in array])


def Ai_tay(array, alpha, beta, n):
    return alpha*arrayize(array, f, n) - beta*arrayize(array, g, n)


def Bi_tay(x, alpha, beta, n):
    return np.sqrt(3)*(alpha*arrayize(x, f, n) + beta*arrayize(x, g, n))


def Ai_pos(x, n):
    x = mpf(x)
    return exp(-xi(x))/(2*sqrt(pi) * x**(mpf(1/4))) * L(-xi(x), n)


def Bi_pos(x, n):
    x = mpf(x)
    return exp(xi(x))/(sqrt(pi) * x**(mpf(1/4))) * L(xi(x), n)


def Ai_neg(x, n):
    x = mpf(x)
    return 1/(sqrt(pi)*(-x)**(1/4)) * (sin(xi(x) - mpf(pi/4))*Q(xi(x), n) + cos(xi(x) - mpf(pi/4))*P(xi(x), n))


def Bi_neg(x, n):
    x = mpf(x)
    return 1/(sqrt(pi)*(-x)**(1/4)) * (-sin(xi(x) - np.pi/4)*P(xi(x), n) + cos(xi(x) - pi/4)*Q(xi(x), n))


a = 0.355028053887817239
b = 0.258819403792806798

x_lin1 = np.linspace(-40, -20, 201)
x_lin2 = np.linspace(-20, 8, 201)  # Taylor je lahko do 8 za Ai, da bo stimala natancnost?
x_lin3 = np.linspace(8, 15, 201)

# plt.plot(x_lin1, arrayize(x_lin1, airyai), c="#ffc07f", label="Ai(x)")
# plt.plot(x_lin1, arrayize(x_lin1, airybi), c="#dac4f7", label="Bi(x)", ls="--")
# plt.plot(x_lin2, arrayize(x_lin2, airyai), c="#ffc07f")
# plt.plot(x_lin2, arrayize(x_lin2, airybi), c="#dac4f7", ls="--")
# plt.plot(x_lin3, arrayize(x_lin3, airyai), c="#ffc07f")
# plt.plot(x_lin3, arrayize(x_lin3, airybi), c="#dac4f7", ls="--")

# plt.plot(x_lin1, arrayize(x_lin1, Ai_neg, 5), c="#8ed081", label=r"$Ai(x)$")
# plt.plot(x_lin1, arrayize(x_lin1, Bi_neg, 5), c="#ffe66d", label=r"$Bi(x)$")
# plt.plot(x_lin2, Ai_tay(x_lin2, a, b, 1000), c="#8ed081")
# plt.plot(x_lin2, Bi_tay(x_lin2, a, b, 1000), c="#ffe66d")
# plt.plot(x_lin3, arrayize(x_lin3, Ai_pos, 10), c="#8ed081")
# plt.plot(x_lin3, arrayize(x_lin3, Bi_pos, 50), c="#ffe66d")
# plt.title("Airyjevi funkciji")
# plt.axhline(alpha=1, ls=":", c="#adadad")
# plt.ylim(-0.6, 0.6)
# plt.legend()
# plt.show()
# plt.show()

# Abs error plot
# x_lin1 = np.linspace(-40, -20, 201)
# x_lin2 = np.linspace(-20, 8, 201)  # Taylor je lahko do 8 za Ai, da bo stimala natancnost?
# x_lin3 = np.linspace(8, 15, 201)
# plt.plot(x_lin1, arrayize(x_lin1, airyai) - arrayize(x_lin1, Ai_neg, 3), c="#540D6E", label="Ai_neg")
# plt.plot(x_lin2, arrayize(x_lin2, airyai) - Ai_tay(x_lin2, a, b, 1000), c="#ee4266", label="Ai_tay")
# plt.plot(x_lin3, arrayize(x_lin3, airyai) - arrayize(x_lin3, Ai_pos, 3), c="#FFD23F", label="Ai_pos")


# plt.plot(x_lin1, arrayize(x_lin1, airybi) - arrayize(x_lin1, Bi_neg, 3), c="#CABAC8", label="Bi_neg")
# plt.plot(x_lin2, arrayize(x_lin2, airybi) - Bi_tay(x_lin2, a, b, 1000), c="#94fbab", label="Bi_tay")
# plt.plot(x_lin3, arrayize(x_lin3, airybi) - arrayize(x_lin3, Bi_pos, 30), c="#f7e733", label="Bi_pos")
#
# plt.title("Absolutna napaka Arijevih funkcij")
# plt.axhline(alpha=1, ls=":", c="#adadad")
# #plt.ylim(-0.6, 0.6)
# plt.yscale("log")
# plt.legend()
# plt.show()

# Rel error plot
x_lin1 = np.linspace(-40, -20, 201)
x_lin2 = np.linspace(-20, 4.5, 201)  # Taylor je lahko do 8 za Ai, da bo stimala natancnost?
x_lin3 = np.linspace(4.5, 15, 201)
# x_lin3 = np.linspace(10**10, 10**11, 201)
plt.plot(x_lin1, (arrayize(x_lin1, airyai) - arrayize(x_lin1, Ai_neg, 3)) / arrayize(x_lin1, airyai), c="#540D6E", label="Ai_neg")
# plt.plot(x_lin2, (arrayize(x_lin2, airyai) - Ai_tay(x_lin2, a, b, 1000)) / arrayize(x_lin2, airyai), c="#ee4266", label="Ai_tay")
# plt.plot(x_lin3, (arrayize(x_lin3, airyai) - arrayize(x_lin3, Ai_pos, 4)) / arrayize(x_lin3, airyai), c="#FFD23F", label="Ai_pos")

# plt.plot(x_lin1, (arrayize(x_lin1, airybi) - arrayize(x_lin1, Bi_neg, 3)) / arrayize(x_lin1, airybi), c="#CABAC8", label="Bi_neg")
# plt.plot(x_lin2, (arrayize(x_lin2, airybi) - Bi_tay(x_lin2, a, b, 1000)) / arrayize(x_lin2, airybi), c="#94fbab", label="Bi_tay")
#plt.plot(x_lin3, (arrayize(x_lin3, airybi) - arrayize(x_lin3, Bi_pos, 300)) / arrayize(x_lin3, airybi), c="#f7e733", label="Bi_pos")

plt.title("Relativna napaka Arijevih funkcij")
plt.axhline(alpha=1, ls=":", c="#adadad")
#plt.ylim(-0.6, 0.6)
plt.yscale("log")
plt.legend()
plt.show()