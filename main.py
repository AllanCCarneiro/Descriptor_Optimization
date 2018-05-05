import numpy as np
from sa import sa
import matplotlib.pyplot as plt
import time
from nmbe_multi import nmbe
from functools import partial
from cost_functions import db, ch, silhouette
import sys

if __name__ == '__main__':
        
    plt.rcParams['text.latex.preamble']=[
    		r"\usepackage{times}",
    		r"\usepackage[brazilian]{babel}",
    		r"\usepackage[utf8x]{inputenc}",
            r"\usepackage[T1]{fontenc}",]
    #Options
    params = {'text.usetex' : True,
              'font.size' : 20,
              'legend.fontsize': 20,
              'font.family' : 'times',
              'text.latex.unicode': True,
    		  'mathtext.default' : 'regular'
              }
    plt.rcParams.update(params)
    
    # Parameters
    
    N_iterations=100 #Number of iterations
    T0=100
    alpha=0.9
    step=0.9
    L=10
    
    im_database_dir='C:/Users/PC/Descriptor_Optimization/pills_balanceada/'
    result_file='C:/Users/PC/Desktop/Descriptor_Optimization/pills_balanceada'
    cost_function=db
    
    search_range=np.array([[0.4,100],[0.4,100],[0.4,100],[0.4,100],[0.4,100]]) #Search space
    
    f=np.zeros([N_iterations])
    
    # The classes.txt file is composed of a python dictionary generated with the pickle module in which the keys are
    # the name of each image file of the dataset and the values corresponding to an array, whose first element is the
    # class that the image belongs to and the others elements correspond to the values of the descriptor for different scales
    descriptor_function=partial(nmbe,database_dir=im_database_dir,database_classes=np.load(im_database_dir+'classes.txt'))
    fitness_function=partial(cost_function,descriptor_function=descriptor_function)
    
    N_int=30
    
    t=np.zeros(N_int)
    f=np.zeros(N_iterations)
    
    
    ax1 = plt.subplot(111)
    
    for i in range(N_int):
        sys.stdout = open(result_file+'_'+cost_function.__name__+'_'+str(i)+'.txt', 'w')
        start = time.time()
        [s,fit_s_vector]=sa(fitness_function,search_range,N_iterations,T0,alpha,L,step)
        f=f+fit_s_vector[1:]
        
        plt.figure(1)
        line1,=plt.plot(-fit_s_vector)
        
        end = time.time()
        
        t[i]=end-start
        
        sys.stdout.close()

    print(end - start)