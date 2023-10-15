#!/usr/bin/env python
import numpy as np
import sys,math

# val = 1 # dot makes all the difference! ... or float()
# print("Variable val : type {} and value {}".format(type(val),val))

# sys.exit(0)

x=float(sys.argv[1])
print ("Vrednost x = {}".format(x))

niter=100
epsilon=1e-6
result=1.
val=1.

for n in range(1,niter+1):
    print ("Iteracija {}".format(n))
    val = val*x/n
    old = result
    result += val
    print ("Korak/clen {:.20e}".format(val))
    print ("Rezultat: stari {:.20f} , novi {:.20f} in korak {:.20f}".format(result,old,val))
    # koncaj prej: 
    # razlicne opcije
    #if val < epsilon:
    #if abs(result-old) == 0.:
    if result == old:
        print("Koncam z rezultatom {:.20f}, tocno {:.20f}, razlika {:.20f}".format(result,math.exp(x),result-math.exp(x)))
        break

