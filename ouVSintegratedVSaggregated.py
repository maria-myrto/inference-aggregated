import numpy as np
import math
import matplotlib.pyplot as plt
import random
random.seed(10)#1000#998
np.random.seed(10)#1000#998

dt = 0.0001

def euler_maruyama(a,sigma,x0,t):
    
    a = float(a)
    sigma = float(sigma)
    print 'drift: '+str(a) 
    print 'diffusion: '+str(sigma) 

    l = t.size
    w = np.zeros(l)
    y = np.zeros(l)
    w[0] = x0
    y[0] = 0 #x0
        
    sqrt_dt = math.sqrt(dt)
    for i in range(0,l-1):
        #### euler-maruyama ##################
        w[i+1] = w[i] -(a)*w[i]*dt + sigma*random.gauss(0,1)*sqrt_dt
        #w[i+1] = w[i] +(a)*(2-w[i])*dt + sigma*random.gauss(0,1)*sqrt_dt
        
        y[i+1] = y[i] + w[i]*dt
    return w , y

def aggregate_measurements(obs_all,dt):
    #print len(obs_all)
    y_aggr = np.zeros(len(obs_all))
    y_aggr[0] = 0#obs_all[0]
    for i in range(0,len(obs_all)-1):
        y_aggr[i+1] = y_aggr[i] + obs_all[i]*dt
    y_aggr_list = y_aggr.tolist()
    return y_aggr, y_aggr_list

n = 20
x0 = 20
step = 2.0
t = np.arange(0,n+dt,dt)
[ou, aggr] = euler_maruyama(4.0,2.0,x0,t)
   
real_obs = x0
#aggr_obs_anal = x0
aggr2 = 0
obs = x0
previous = 0
list_aggr_all = []
#obs2 = x0
for  i in range(0,int((1/step)*n)):
    obs_all = ou[previous:(float(i+1)*step)*(1/dt)+1]  
    [y_aggr_2, y_list] = aggregate_measurements(obs_all,dt)
    list_aggr_all.append(y_list)
    aggr2 = np.vstack((aggr2,y_aggr_2[len(y_aggr_2)-1]))         
    real_obs = np.vstack((real_obs,ou[(float(i+1)*step)*(1/dt)]))
    aggr_obs = dt*0.5*(ou[(float(i)*step)*(1/dt)]+2*(sum(ou[((float(i)*step)*(1/dt)+1):((float(i+1)*step)*(1/dt))]))+ou[((float(i+1)*step)*(1/dt))])
    obs = np.vstack((obs,aggr_obs))
    previous = (float(i+1)*step)*(1/dt)+1   
aggr_all = np.concatenate(list_aggr_all)
   



plt.figure(1)
plt.plot(t,ou)
plt.xlabel('time',fontsize=16)
plt.ylabel('$X_t$',fontsize=16)
#plt.title('OU',fontsize=16)

plt.figure(2)
plt.plot(t,aggr)    
plt.xlabel('time',fontsize=16)
plt.ylabel('$Y_t$',fontsize=16)
#plt.title('Integral of OU',fontsize=16)
    
plt.figure(3)
plt.plot(t,aggr_all)    
plt.xlabel('time',fontsize=16)
plt.plot(np.arange(0,n+step,step),obs,'rx')
plt.ylabel('$Y_t$',fontsize=16)
plt.ylim((-2,5))
#plt.title('Aggregated process',fontsize=16)


