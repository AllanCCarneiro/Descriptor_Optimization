# -*- coding: utf-8 -*-
# descriptors: module that implements the calculation of signatures and descriptors of images

import numpy as np
import cv2 as cv
from scipy.interpolate import interp1d 
from scipy.spatial.distance import pdist,squareform
from math import sqrt,acos

class contour_base:
 '''Represents an binary image contour as a complex discrete signal. 
   Some calculation methods are provided to compute contour 1st derivative, 2nd derivatives and perimeter.
   The signal variable (self.c) is represented as a single dimensional ndarray of complex.  
   This class is interable and callable so, interation over objects results in sequential access to each signal variable element. Furthermore, calling the object as a function yields as return value te signal variable.
 
 '''

 def __init__(self,fn):

  self.__i = 0
  if type(fn) is str:
   im = cv.imread(fn,cv.IMREAD_GRAYSCALE)
   im = cv.imread(fn)
   imgray = cv.cvtColor(im, cv.COLOR_BGR2GRAY)
   ret, thresh = cv.threshold(imgray, 127, 255, 0)
   im2, s, hierarchy = cv.findContours(thresh,cv.RETR_LIST,cv.CHAIN_APPROX_NONE)
   s=np.squeeze(s[0])
   self.c = np.array([complex(i[1],i[0]) for i in s])
  elif (type(fn) is np.ndarray):
    self.c = fn
  elif (type(fn) is cv.iplimage):
    s = cv.findContours(fn,cv.RETR_LIST,cv.CHAIN_APPROX_NONE) 
    self.c = np.array([complex(i[1],i[0]) for i in s])
  N = self.c.size
  self.freq = np.fft.fftfreq(N,1./float(N))

  self.ftc = np.fft.fft(self.c)

  if isinstance(self,contour_base):
   self.calc_derivatives()
 
 def calc_derivatives(self):
   ftcd = np.complex(0,1) * 2 * np.pi * self.freq * self.ftc
   ftcdd = - (2 * np.pi * self.freq)**2 * self.ftc 
   self.cd = np.fft.ifft(ftcd)
   self.cdd = np.fft.ifft(ftcdd)

 def first_deriv(self): 
  '''Return the contour signal 1st derivative'''
  return self.cd

 def second_deriv(self): 
  '''Return the contour signal 2nd derivative'''
  return self.cdd

 def perimeter(self): 
  '''Calculate and return the contour perimeter'''
  return (2*np.pi*np.sum(np.abs(self.cd))/float(self.cd.size))

 def __iter__(self): return self

 def next(self):

  if self.__i > self.c.size-1:
   self.__i = 0
   raise StopIteration
  else:
   self.__i += 1
   return self.c[self.__i-1]
 
 def __call__(self): return self.c

class contour(contour_base):
  '''Like contour_base except that, prior to derive a complex signal representation, smooths the image contour using a Gaussian kernel. The kernel parameter (gaussian standard deviation) is the second constructor parameter. See also contour_base.'''

  # Gaussian smoothing function 
  def __G(self,s):
    return (1/(s*(2*np.pi)**0.5))*np.exp(-self.freq**2/(2*s**2))
  
  def __init__(self,fn,sigma=None):
   contour_base.__init__(self,fn)
   if sigma is not None:
    E = np.sum(self.ftc * self.ftc.conjugate())
    self.ftc = self.ftc * self.__G(sigma)
    Eg  = np.sum(self.ftc * self.ftc.conjugate())
    k = sqrt(abs(E/Eg))
    self.c = np.fft.ifft(self.ftc)*k
    self.calc_derivatives()
    self.cd = self.cd * k
    self.cdd = self.cdd * k
     

# class curvature: calculates the curvature of a contour for various levels of smoothing
# constructor parameters:   def __init__(self,fn = None,sigma_range = np.linspace(2,30,10)) 
# fn : It can be the name of an image file (string) that contains a binary form or a vector (ndarray) of values of the
# contour coordinates of a shape (complex representation x + j.y).
# In the first case, the contours are extracted through the cv.FindContours() function of the Opencv library
# sigma_range: vector (ndarray) containing the values that will be used as the standard deviation for the Gaussian PDF.
# which filters the contours before calculating the curvature.
#  when zero no filtering is applied to contour

class curvatura:
  '''For a given binary image calculates and yields a family of curvature signals represented in a two dimensional ndarray structure; each row corresponds to the curvature signal derived from the smoothed contour for a certain smooth level.'''

  def __Calcula_Curvograma(self,fn):
   z = contour(fn)
   caux = [contour(z(),s) for s in self.sigmas]
   caux.append(z)
   self.contours = np.array(caux)
   self.t = np.linspace(0,1,z().size)
   self.curvs = np.ndarray((self.sigmas.size+1,self.t.size),dtype = "float")
  
   for c,i in zip(self.contours,np.arange(self.contours.size)):
    # Calculates curvature for various smoothing scales of shape contour
     curv = c.first_deriv() * np.conjugate(c.second_deriv())
     curv = - curv.imag
     curv = curv/(np.abs(c.first_deriv())**3)
     # two-dimensional array curvs = Curvature Function k(sigma,t) 
     self.curvs[i] = np.copy(np.tanh(curv))   
 
  # Contructor
  def __init__(self,fn = None,sigma_range = np.linspace(2,30,20)):
   # Extract image contour
   self.sigmas = sigma_range
   self.__Calcula_Curvograma(fn)

 # Function to compute curvature
 # It is called into class constructor
  def __call__(self,idx = 0,t= None):
    if t is None:
     __curv = self.curvs[idx]
    elif (type(t) is np.ndarray):
     __curv = interp1d(self.t,y = self.curvs[idx],kind='quadratic')
     return(__curv(t))
    else:
      __curv = self.curvs[idx]

    return(__curv)
   
class bendenergy:
 ''' For a given binary image, computes the multiscale contour curvature bend energy descriptor'''
 
 def __init__(self,fn,scale):
  self.__i = 0
  k = curvatura(fn,scale[::-1])
  # p = non-smoothed contour perimeter
  p = k.contours[-1].perimeter() 
  self.phi  = np.array([(p**2)*np.mean(k(i)**2) for i in np.arange(0,scale.size)])

 def __call__(self): return self.phi
 
 def __iter__(self): return self

 def next(self):

   if self.__i > self.phi.size-1:
    self.__i = 0
    raise StopIteration
   else:
    self.__i += 1
    return self.phi[self.__i-1]


