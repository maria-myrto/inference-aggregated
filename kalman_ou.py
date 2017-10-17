import numpy as np
import math
import scipy
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from ou_odes import ou_odes
from ou_simulation_1D import ou_simulation_1D
import random 

random.seed(100)
np.random.seed(100)

n = 10
x0 = 20
step = 1.0

real_theta = 4.0
real_sigma = 2.0

[obs_trapz, ou, time1, obs_anal,trace] = ou_simulation_1D(x0, n, step,real_theta,real_sigma)
real_obs = ou[1:len(ou)]
obs = obs_anal[1:len(obs_anal)]
#obs = real_obs
#noise = np.matrix(np.random.normal(0, 0.01, len(obs)))  # .T
obs = obs#+noise.T
obs_avg = obs/step
#obs_avg = real_obs
#plt.figure(0)
#plt.plot(obs_anal,'ko')
#plt.plot(obs_trapz,'ro')
#plt.plot(ou,'go')

#theta = 4.0
#sigma = 2.0

c1_mat = np.load("C:\data_ou\c1_mat1_all1.dat")
c1_post_mean = np.mean(c1_mat)
theta = np.mean(c1_post_mean)

c2_mat = np.load("C:\data_ou\c1_mat1_all1.dat")
c2_post_mean = np.mean(c2_mat)
sigma = np.mean(c2_post_mean)

P = np.matrix([1])
V = np.matrix([0])
m0 = np.matrix([x0])
S0 = np.matrix([0.5])
P0 = np.matrix([1])
V0 = np.matrix([0])

#c = np.arry([theta,mu,sigma])
c = np.array([theta,sigma])

mean0 = P0*m0.T
cov0 = P0*S0*P0.T+V0

y0 = np.random.multivariate_normal(np.array(mean0.T)[0], np.array(cov0))
y0 = np.matrix(y0)
#y0 = np.matrix([x0])
K_0 = S0*P0.T*(P0*S0*P0.T+V0).I
m0_star = m0.T + K_0*(y0.T-P0*m0.T)
S0_star = S0 - K_0*P0*S0

init = np.array([m0_star[0,0],S0_star[0,0]])
y_det = np.array([0])
y_all = np.array([0,0.1])
t_all = np.array([0])

j = 0

dt = 0.0001

m_star_mat = np.array([m0_star[0,0]])
S_star_mat = np.array([S0_star[0,0]])
for i in np.arange(0,n,step):
    time = np.arange(i,i+step+dt,dt)
    t0 = time[0]
    t1 = time[len(time)-1]
    delta = t1-t0
#   ################### numerical solution of odes ################################
    y = odeint(ou_odes,init,time,args = (c,))
    l = len(y)
    t_all = np.concatenate((t_all,time))
#    y_det = np.vstack((y_det,y[:,0:1]))
    y_all = np.vstack((y_all,y))
    m = np.matrix([y[l-1,0]])
    S = np.matrix([y[l-1,1]])
###########################################################    
##    mean = P*m.T
##    cov = P*S*P.T+V
##    det = np.linalg.det(cov)
##    a = (np.power((2*np.pi),-1))*np.power(det,-0.5)
##    b = np.exp(-0.5*(obs[j]-mean)*cov.I*(obs[j]-mean).T)
##    lik = a*b
    
########################## analytical solution of odes c, tau #########################    
#    m_anal = init[0]*np.exp(-(time[len(time)-1]-time[0])/c[0])
#    S_anal = (float(c[0])*c[1]/2)*(1- np.exp(-2*(time[len(time)-1]-time[0])/c[0]))+init[1]*np.exp(-2*(time[len(time)-1]-time[0])/c[0])
#
#    m = np.matrix([[m_anal]])
#    S = np.matrix([[S_anal]])
#    #print S
######################### analytical solution of odes a, lamda #########################    
#    #m_anal = init[0]*np.exp(-(time[len(time)-1]-time[0])*c[0])
#    #S_anal = ((c[1]**2)/(float(c[0])*2))*(1- np.exp(-2*(time[len(time)-1]-time[0])*c[0]))+init[1]*np.exp(-2*(time[len(time)-1]-time[0])*c[0])
#    m_anal = 1.0*init[0]*np.exp(-1.0*delta*c[0])
#    S_anal = ((c[1]**2)/(c[0]*2.0))*(1.0- np.exp(-2.0*delta*c[0]))+init[1]*np.exp(-2.0*delta*c[0])
#  
#    m = np.matrix([[m_anal]])
#    S = np.matrix([[S_anal]])
#
#################################################################################################    
    K = S*P.T*(P*S*P.T+V).I
    m_star = m.T + K*(obs_avg[j].T-P*m.T)
    m_star_mat = np.vstack((m_star_mat,np.array(m_star).T))
    S_star = S - K*P*S
    S_star_mat = np.vstack((S_star_mat,np.array(S_star)))
    
#    print m_star
#    print S_star
#    x[j,:] = np.random.multivariate_normal(m_star,S_star)
    init = np.array([m_star[0,0],S_star[0,0]])
    #init = np.array([m_star[0,0],0])
    j +=1

mean = m_star_mat[1:len(m_star_mat)]

variance = S_star_mat[1:len(S_star_mat)]


#xronos = np.arange(0,n+step,step)
xronos = np.arange(step,n+step,step)
# plt.ion()
plt.figure(5)
plt.plot(time1,trace)
plt.plot(xronos,obs_avg,'ro')
#plt.plot(xronos,obs,'ko')

#plt.plot(xronos,mean1)
#plt.xlim(0,2.5)
#plt.plot(xronos,mean-np.sqrt(variance),'g')
#plt.plot(xronos,mean+np.sqrt(variance),'g')
#plt.plot(xronos,mean)

#plt.xlim(3,10)
#plt.ylim(-5,20)
plt.plot(t_all,y_all[:,0],'k')
plt.plot(t_all,y_all[:,0]-np.sqrt(y_all[:,1]),'g')
plt.plot(t_all,y_all[:,0]+np.sqrt(y_all[:,1]),'g')
plt.title('KF1-OU')
##plt.plot(xronos,obs,'go')
##plt.draw()
plt.show()

#plt.figure(6)
#plt.plot(xronos,real_obs,'r')
#plt.plot(xronos,obs,'g')