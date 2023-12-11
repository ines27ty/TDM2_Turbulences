def langevin():
    # This python function write a txt file with a random signal obeying a langevin stochastic equation
    
    import numpy as np
    dt = 0.01 #the time step of the signal
    Nit = 10000 #the number of time step of the signal 
    
    sigma = 2  #the variance of the langevin signal
    T = 1.  #tha characteristic time of the signal
    mu = 3.  #the mean value of the signal
    

    sigma2 = sigma*sigma  
    k= np.sqrt(2.*dt/T) 

    #the intial condition
    dw = np.random.randn()
    x0 = sigma * np.random.randn()+mu
    t0 = 0.
    
    Y = [t0,x0,dw]
    
    i=0
    x=x0
    t=t0
    while i<Nit:
        i += 1
        dw = k * np.random.randn()
        dx = sigma * dw
        dx += (mu-x)*dt/T
        
        x += dx
        t += dt
        Y = np.vstack( (Y, [t,x,dw]) )  
         
    np.savetxt("langevin.txt",Y) 