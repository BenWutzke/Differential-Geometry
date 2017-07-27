# -*- coding: utf-8 -*-

import math
import scipy.integrate as integrate
import numpy as np
import matplotlib.pyplot as plt

def curve(a,b):
   X = lambda t: sum([a[i]*t**i for i,_ in enumerate(a)])
   Y = lambda t: sum([b[i]*t**i for i,_ in enumerate(b)])
   return X,Y

def curve_prime(a,b):
   X = lambda t: sum([i*a[i]*t**(i-1) if i not in [0] else 0 for i,_ in enumerate(a)])
   Y = lambda t: sum([i*b[i]*t**(i-1) if i not in [0] else 0 for i,_ in enumerate(b)])
   return X,Y

def curve_pprime(a,b):
   X = lambda t: sum([i*(i-1)*a[i]*t**(i-2) if i not in [0,1] else 0 for i,_ in enumerate(a)])
   Y = lambda t: sum([i*(i-1)*b[i]*t**(i-2) if i not in [0,1] else 0 for i,_ in enumerate(b)])
   return X,Y

def curvature(a,b):
   X_p,Y_p = curve_prime(a,b)
   X_pp,Y_pp = curve_pprime(a,b)
   return lambda t: math.pow(X_p(t)**2 + Y_p(t)**2,3/2.0)/math.sqrt(X_p(t)*Y_pp(t) + Y_p(t)*X_pp(t))

def curve_length(a,b,T):
   X_p,Y_p = curve_prime(a,b)
   return lambda t: integrate.quad(math.sqrt(X_p(t) + Y_p(t)),T[1],T[2])

def curve_length_prime(a,b):
   X_p,Y_p = curve_prime(a,b)
   X_pp,Y_pp = curve_pprime(a,b)
   return lambda t: 0 if X_p(t) + Y_p(t) <=0 else (X_pp(t) + Y_pp(t))/math.sqrt(X_p(t) + Y_p(t))

def tangential_angle(a,b):
   K = curvature(a,b)
   s_p = curve_length_prime(a,b)
   ta = lambda t: s_p(t)*K(t)
   return lambda t: integrate.quad(ta,0,t)[0]

def main():
   """
   We define a parametric curve here, where the values of a and b correspond to to the coefficients of the polynomials that describe 
   x and y space.
   x(t) = a0t^0 + a1t^1 + a2t^2 + ...
   y(t) = b0t^0 + b1t^1 + b2t^2 + ...
   """
   a = [1,1,-5]
   b = [1,-5,-1]

   x,y = curve(a,b)
   tan_ang = tangential_angle(a,b)
   T = np.arange(-1,1,0.1)
   plt.plot([x(t) for t in T],[y(t) for t in T])
   #plt.plot(T,[tan_ang(t) for t in T])
   plt.plot()
   plt.show()      
   for t in T:
      print(str(x(t))," ",str(y(t))," ",str(tan_ang(t)))
      
if __name__ == '__main__':
   main()
