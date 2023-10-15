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


def f_rec(x, n):
    x = mpf(x)
    c = rf(1 / 3, 0) * mpf("3") ** mpf(0) * mpf(x) ** mpf(3 * 0) / fac(3 * 0)
    for k in range(0, n+1):
        if k == 0:
            prevstep = c
            continue
        print("f: " + str(k))
        step = prevstep * (mpf(1/3) + k - 1) * 3 * x**3/((3*k)*(3*k - 1)*(3*k - 2))
        old = c
        c += step
        prevstep = step
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


def g_rec(x, n):
    x = mpf(x)
    c = rf(2 / 3, 0) * mpf("3") ** mpf(0) * mpf(x) ** mpf(3 * 0 + 1) / fac(3 * 0 + 1)
    for k in range(0, n+1):
        if k == 0:
            prevstep = c
            continue
        print("g: " + str(k))
        # step = prevstep * rf(2/3, k) * mpf("3")**mpf(k) * mpf(x)**mpf(3*k + 1)/fac(3*k + 1)
        step = prevstep * (mpf(2/3) + k - 1) * 3 * x**3/((3*k + 1)*(3*k)*(3*k - 1))
        old = c
        c += step
        prevstep = step
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
    for s in range(0, n+1):
        print("L:" + str(s))
        step = u_s(s)/x**s
        c += step
        if s == 0:
            prevstep = mpf(10e100)  # Samo nekaj zelo velikega, da ima na zacetku za primerjat z neko neskoncnostjo
            continue
        if step >= prevstep:
            return c
        prevstep = step
    return c


def P(x, n):
    c = mpf(0)
    x = mpf(x)
    for s in range(0, n+1):
        print("P: " + str(s))
        step = (-1)**s * u_s(2*s)/mpf(x)**(2*s)
        c += step
        if s == 0:
            prevstep = mpf(10e100)  # Samo nekaj zelo velikega, da ima na zacetku za primerjat z neko neskoncnostjo
            continue
        if step >= prevstep:
            return c
        prevstep = step
    return c


def Q(x, n):
    c = mpf(0)
    x = mpf(x)
    for s in range(0, n+1):
        print("Q: " + str(s))
        step = (-1)**s * u_s(2*s + 1) / mpf(x) ** (2*s + 1)
        c += step
        if s == 0:
            prevstep = mpf(10e100)
            continue
        if step >= prevstep:
            return c
    return c


def arrayize(array, function, *args):
    return np.array([function(x, *args) for x in array])


def Ai_tay(array, alpha, beta, n):
    return alpha*arrayize(array, f, n) - beta*arrayize(array, g, n)


def Ai_tay_rec(array, alpha, beta, n):
    return alpha*arrayize(array, f_rec, n) - beta*arrayize(array, g_rec, n)


def Bi_tay(x, alpha, beta, n):
    return np.sqrt(3)*(alpha*arrayize(x, f, n) + beta*arrayize(x, g, n))


def Bi_tay_rec(x, alpha, beta, n):
    return np.sqrt(3)*(alpha*arrayize(x, f_rec, n) + beta*arrayize(x, g_rec, n))


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
#
# plt.title("Airyjevi funkciji")
# plt.axhline(alpha=1, ls=":", c="#adadad")
# plt.ylim(-0.6, 0.6)
# plt.legend()
# plt.xlabel("x")
# plt.ylabel(r"$f(x)$")
# plt.show()

# plt.plot(x_lin1, arrayize(x_lin1, Ai_neg, 5), c="#8ed081", label=r"$Ai(x)$")
# plt.plot(x_lin1, arrayize(x_lin1, Bi_neg, 5), c="#ffe66d", label=r"$Bi(x)$")
# plt.plot(x_lin2, Ai_tay(x_lin2, a, b, 1000), c="#8ed081")
# #plt.plot(x_lin2, Ai_tay_rec(x_lin2, a, b, 1000), c="#f76f8e")
# plt.plot(x_lin2, Bi_tay(x_lin2, a, b, 1000), c="#ffe66d")
# #plt.plot(x_lin2, Bi_tay_rec(x_lin2, a, b, 1000), c="#ff69eb")
# plt.plot(x_lin3, arrayize(x_lin3, Ai_pos, 10), c="#8ed081")
# plt.plot(x_lin3, arrayize(x_lin3, Bi_pos, 50), c="#ffe66d")
# plt.title("Airyjevi funkciji")
# plt.axhline(alpha=1, ls=":", c="#adadad")
# plt.ylim(-0.6, 0.6)
# plt.legend()
# plt.xlabel("x")
# plt.ylabel(r"$f(x)$")
# plt.show()

# Abs error plot
# x_lin1 = np.linspace(-40, -20, 201)
# x_lin2 = np.linspace(-20, 8, 201)  # Taylor je lahko do 8 za Ai, da bo stimala natancnost?
# x_lin3 = np.linspace(8, 15, 201)
plt.plot(x_lin1, np.abs(arrayize(x_lin1, airyai) - arrayize(x_lin1, Ai_neg, 3)), c="#540D6E", label=r"$Ai_{neg}$")
plt.plot(x_lin2, np.abs(arrayize(x_lin2, airyai) - Ai_tay(x_lin2, a, b, 1000)), c="#ee4266", label=r"$Ai_{tay}$")
plt.plot(x_lin3, np.abs(arrayize(x_lin3, airyai) - arrayize(x_lin3, Ai_pos, 30)), c="#FFD23F", label=r"$Ai_{pos}$")

plt.title("Absolutna napaka Arijevih funkcij")
plt.axhline(alpha=1, ls=":", c="#adadad")
plt.yscale("log")
plt.legend()
plt.xlabel("x")
plt.ylabel(r"$|Ai_{ref} - f|$")
plt.show()

plt.plot(x_lin1, np.abs(arrayize(x_lin1, airybi) - arrayize(x_lin1, Bi_neg, 30)), c="#CABAC8", label=r"$Bi_{neg}$")
plt.plot(x_lin2, np.abs(arrayize(x_lin2, airybi) - Bi_tay(x_lin2, a, b, 1000)), c="#94fbab", label=r"$Bi_{tay}$")
plt.plot(x_lin3, np.abs(arrayize(x_lin3, airybi) - arrayize(x_lin3, Bi_pos, 30)), c="#f7e733", label=r"$Bi_{pos}$")

plt.title("Absolutna napaka Arijevih funkcij")
plt.axhline(alpha=1, ls=":", c="#adadad")
plt.yscale("log")
plt.legend()
plt.xlabel("x")
plt.ylabel(r"$|Bi_{ref} - f|$")
plt.show()

# Rel error plot
#x_lin1 = np.linspace(-30, -20, 201)
#x_lin2 = np.linspace(-20, 4.5, 201)  # Taylor je lahko do 8 za Ai, da bo stimala natancnost?
#x_lin3 = np.linspace(4.5, 15, 201)
# x_lin3 = np.linspace(10**10, 10**11, 201)
plt.plot(x_lin1, np.abs((arrayize(x_lin1, airyai) - arrayize(x_lin1, Ai_neg, 30))) / arrayize(x_lin1, airyai), c="#540D6E", label=r"$Ai_{neg}$")
plt.plot(x_lin2, np.abs((arrayize(x_lin2, airyai) - Ai_tay(x_lin2, a, b, 1000))) / arrayize(x_lin2, airyai), c="#ee4266", label=r"$Ai_{tay}$")
plt.plot(x_lin3, np.abs((arrayize(x_lin3, airyai) - arrayize(x_lin3, Ai_pos, 4))) / arrayize(x_lin3, airyai), c="#FFD23F", label=r"$Ai_{pos}$")

plt.title("Relativna napaka Arijevih funkcij")
plt.axhline(alpha=1, ls=":", c="#adadad")
plt.yscale("log")
plt.legend()
plt.ylabel(r"$\frac{|Ai_{ref} - f|}{Ai_{ref}}$")
plt.xlabel("x")
plt.show()

plt.plot(x_lin1, np.abs((arrayize(x_lin1, airybi) - arrayize(x_lin1, Bi_neg, 3)) / arrayize(x_lin1, airybi)), c="#CABAC8", label=r"$Bi_{neg}$")
plt.plot(x_lin2, np.abs((arrayize(x_lin2, airybi) - Bi_tay(x_lin2, a, b, 1000))) / arrayize(x_lin2, airybi), c="#94fbab", label=r"$Bi_{tay}$")
plt.plot(x_lin3, np.abs((arrayize(x_lin3, airybi) - arrayize(x_lin3, Bi_pos, 300))) / arrayize(x_lin3, airybi), c="#f7e733", label=r"$Bi_{pos}$")

plt.title("Relativna napaka Arijevih funkcij")
plt.axhline(alpha=1, ls=":", c="#adadad")
plt.yscale("log")
plt.legend()
plt.ylabel(r"$\frac{|Bi_{ref} - f|}{Bi_{ref}}$")
plt.xlabel("x")
plt.show()