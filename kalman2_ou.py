import numpy as np
import math
import scipy
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from ou_odes2 import ou_odes2
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


plt.figure(0)
plt.plot(obs_anal,'ko')
plt.plot(obs_trapz,'ro')
plt.plot(ou,'go')

#obs = np.load("observations.dat")
#real_obs = np.load("real_observations.dat")

P = np.matrix([1])
V = np.matrix([0])
m0 = np.matrix([x0])
S0 = np.matrix([0.5])
P0 = np.matrix([1])
V0 = np.matrix([0])


#theta = 4.0
#sigma = 2.0

c1_mat = np.load("C:\data_ou\c1_mat2_all01.dat")
c1_post_mean = np.mean(c1_mat)
theta = np.mean(c1_post_mean)

c2_mat = np.load("C:\data_ou\c2_mat2_all01.dat")
c2_post_mean = np.mean(c2_mat)
sigma = np.mean(c2_post_mean)


#mu = 0
#c = np.array([theta,mu,sigma])
c = np.array([theta,sigma])
##mean0 = P0*m0.T
##cov0 = P0*S0*P0.T+V0
##y0 = np.random.multivariate_normal(np.array(mean0.T)[0],np.array(cov0))
y0 = np.matrix(x0)
##det0 = np.linalg.det(cov0)
##a0 = (np.power((2*np.pi),-1))*np.power(det0,-0.5)
##b0 = np.exp(-0.5*(y0-mean0)*cov0.I*(y0-mean0).T)
##lik0 = a0*b0

K_0 = S0*P0.T/(P0*S0*P0.T+V0)
m0_star = m0.T + K_0*(y0.T-P0*m0.T)
S0_star = S0 - K_0*P0*S0
##print y0
##print lik0
##print m0_star
##print S0_star

#init = np.array([m0_star[(0,0)],m0_star[(1,0)],S0_star[(0,0)],S0_star[(0,1)],S0_star[(1,1)],m0_star[(0,0)],m0_star[(1,0)],0,0,0,0,S0_star[(0,0)],S0_star[(0,1)],S0_star[(1,1)]])
init = np.array([m0_star[(0,0)],S0_star[(0,0)],0*m0_star[(0,0)],0*S0_star[(0,0)],0*S0_star[(0,0)]])
#init = np.array([m0_star[(0,0)],0,m0_star[(0,0)],0,2*0])

y_all = np.array([0,0.1,0,0.1])
t_all = np.array([0])

j = 0
#step = 0.5
#x = np.zeros((5,2))
dt = 0.0001
mx_star_mat = np.array([m0_star[0,0]])
S_star_mat = np.array([S0_star[0,0]])
m_star_mat = np.array([m0_star[0,0]])
Q_star_mat = np.array([S0_star[0,0]])
for i in np.arange(0,n,step):
    
    time = np.arange(i, i+step+dt, dt)
    t0 = time[0]
    t1 = time[len(time)-1]
    delta = t1-t0
#    ############################# analytical solution of odes ###############################
#    mx_anal = init[0]*np.exp(-(time[len(time)-1]-time[0])/c[0])
#    ###S_anal = (float(c[0])*c[1]/2)*(1- np.exp(-2*(time[len(time)-1]-time[0])/c[0]))
#    S_anal = (float(c[0])*c[1]/2)*(1- np.exp(-2*(time[len(time)-1]-time[0])/c[0]))+init[1]*np.exp(-2*(time[len(time)-1]-time[0])/c[0])
#    m_anal = init[2]+init[0]*c[0]*(1-np.exp(-(time[len(time)-1]-time[0])/(c[0])))
#        ###QM_anal = (float(c[1])*c[0]**2/2) + (init[3]+init[0]**2*c[0]-c[1]*c[0]**2)*np.exp(-(time[len(time)-1]-time[0])/(c[0])) + (float(c[1])*c[0]**2/2-init[0]**2*c[0])*np.exp(-2*(time[len(time)-1]-time[0])/(c[0]))   
##    ##QM_anal = c[1]*(c[0]**2)*0.5*(1-2*np.exp(-(time[len(time)-1]-time[0])/(c[0]))+np.exp(-2*(time[len(time)-1]-time[0])/c[0]))        
#    QM_anal = c[1]*(c[0]**2)*0.5+(-c[1]*(c[0]**2)+init[1]*c[0]+0*init[1])*np.exp(-(time[len(time)-1]-time[0])/(c[0]))+(c[1]*(c[0]**2)*0.5-init[1]*c[0])*np.exp(-2*(time[len(time)-1]-time[0])/c[0])
#    ###QM_anal = c[1]*(c[0]**2)*0.5+(-c[1]*(c[0]**2)+init[1]*c[0]+init[1]-init[0]*init[2])*np.exp(-(time[len(time)-1]-time[0])/(c[0]))+(c[1]*(c[0]**2)*0.5-init[1]*c[0])*np.exp(-2*(time[len(time)-1]-time[0])/c[0])    
##    ##Q_anal = c[1]*c[0]**3*((time[len(time)-1]-time[0])/(c[0])-2*(1- np.exp(-(time[len(time)-1]-time[0])/(c[0])))+0.5*(1- np.exp(-2*(time[len(time)-1]-time[0])/(c[0]))))
#    Q_anal = c[1]*c[0]**2*(time[len(time)-1]-time[0])+0*init[1]-2*(c[1]*c[0]**3-init[1]*c[0]**2-0*init[1]*c[0])*(1- np.exp(-(time[len(time)-1]-time[0])/(c[0])))+(0.5*c[1]*c[0]**3-init[1]*c[0]**2)*(1- np.exp(-2*(time[len(time)-1]-time[0])/(c[0])))
#    ###Q_anal = init[1]**2-init[0]*init[2]+c[1]*c[0]**2*(time[len(time)-1]-time[0])/(c[0])-2*(c[1]*c[0]**3*+init[1]*c[0]**2)*(1- np.exp(-(time[len(time)-1]-time[0])/(c[0])))+(0.5*c[1]*c[0]**3*-init[1]*c[0]**2)*(1- np.exp(-2*(time[len(time)-1]-time[0])/(c[0])))
############################################ analytical a , s##################################
#    mx_anal = 1.0*init[0]*np.exp(-1.0*delta*c[0])
#    S_anal = ((c[1]**2)/(c[0]*2.0))*(1.0- np.exp(-2.0*delta*c[0]))+init[1]*np.exp(-2.0*delta*c[0])
#        
#    m_anal = 1.0*init[0]/c[0]*(1.0-np.exp(-1.0*delta*c[0]))        
#        #QM_anal = c[1]*(c[0]**2)*0.5+(-c[1]*(c[0]**2)+init[1]*c[0]+init[1])*np.exp(-(time[len(time)-1]-time[0])/(c[0]))+(c[1]*(c[0]**2)*0.5-init[1]*c[0])*np.exp(-2*(time[len(time)-1]-time[0])/c[0])
#        #QM_anal = 0.5*((c[1]**2.0)/(c[0]**2))+(-(c[1]**2.0)/(c[0]**2.0)+1.0*init[1]/c[0])*np.exp(-1.0*delta*(c[0]))+(0.5*((c[1]**2.0)/(c[0]**2.0))-1.0*init[1]/c[0])*np.exp(-2.0*delta*c[0])
#    QM_anal = 0.5*((c[1]**2.0)/(c[0]**2.0))+(-(c[1]**2.0)/(c[0]**2.0)+1.0*init[1]/c[0])*np.exp(-1.0*delta*(c[0]))+(0.5*((c[1]**2.0)/(c[0]**2.0))-1.0*init[1]/c[0])*np.exp(-2.0*delta*c[0])
#
#        #Q_anal = c[1]*(c[0]**2)*(time[len(time)-1]-time[0])+init[1]-2*(c[1]*c[0]**3-init[1]*c[0]**2-init[1]*c[0])*(1- np.exp(-(time[len(time)-1]-time[0])/(c[0])))+(0.5*c[1]*c[0]**3-init[1]*c[0]**2)*(1- np.exp(-2*(time[len(time)-1]-time[0])/(c[0])))
#        #Q_anal = ((c[1]**2.0)/((float(c[0])**2)))*delta*1.0-2.0*((c[1]**2.0)/c[0]**3.0-1.0*init[1]/(c[0]**2.0))*(1.0- np.exp(-delta*c[0]*1.0))+(0.5*(c[1]**2.0)/(c[0]**3.0)-init[1]*1.0/(c[0]**2.0))*(1.0- np.exp(-2.0*delta*(c[0])))
#    Q_anal = 1.0*((c[1]**2.0)/((float(c[0])**2)))*delta-2.0*((c[1]**2.0)/c[0]**3.0-1.0*init[1]/(c[0]**2.0))*(1.0- np.exp(-1.0*delta*c[0]))+(0.5*(c[1]**2.0)/(c[0]**3.0)-1.0*init[1]/(c[0]**2.0))*(1.0- np.exp(-2.0*delta*(c[0])))
#
######################################################
#    mx = np.matrix([[mx_anal]])
#    S = np.matrix([[S_anal]])
#    m = np.matrix([[m_anal]])
#    QM = np.matrix([[QM_anal]])
#    Q = np.matrix([[Q_anal]])
#    #print S_anal, Q_anal, QM_anal  

###################### numerical solution of odes ############################################
    y = odeint(ou_odes2,init,time,args = (c,))
    l = len(y)
    t_all = np.concatenate((t_all,time))
    y_all = np.vstack((y_all,y[:,0:4]))
    mx = np.matrix([y[l-1,0]])    
    S = np.matrix([y[l-1,1]])    
    m = np.matrix([y[l-1,2]])    
    Q2 = np.matrix([y[l-1,4]])
    Q = Q2-m**2    
    mQM = np.matrix([y[l-1,3]])
    QM = mQM-m*mx
###########################################################    
##    mean = P*m.T
##    cov = P*Q*P.T+V
##    det = np.linalg.det(cov)
##    a = (np.power((2*np.pi),-1))*np.power(det,-0.5)
##    b = np.exp(-0.5*(obs[j]-mean)*cov.I*(obs[j]-mean).T)
##    lik = a*b

    part1 = Q*P.T/(P*Q*P.T+V)
    m_star = m.T + part1*(obs[j].T-P*m.T)
    m_star_mat = np.vstack((m_star_mat,np.array(m_star).T))
    
    Q_star = Q - part1*P*Q
    Q_star_mat = np.vstack((Q_star_mat,np.array(Q_star)))

    part1x = QM.T*P.T/(P*Q*P.T+V)
    #print part1x
    mx_star = mx.T + part1x*(obs[j].T-P*m.T)
    #part1x = QM_anal.T*P.T*(P*Q_anal*P.T+V).I
    #mx_star = mx_anal.T + part1x*(obs[j].T-P*m_anal.T)
    mx_star_mat = np.vstack((mx_star_mat,np.array(mx_star).T))
    
    S_star = S - part1x*P*QM
    #S_star = S_anal - part1x*P*QM_anal 
    S_star_mat = np.vstack((S_star_mat,np.array(S_star)))
    
    #print mx_star
    #print S_star
    #x[j,:] = np.random.multivariate_normal([mx_star[0,0],mx_star[1,0]],S_star)
    
    init = np.array([mx_star[(0,0)],S_star,0*mx_star[(0,0)],0*S_star,0*S_star])
    #init = np.array([mx_star[(0,0)],S_star,mx_star[(0,0)],S_star+mx_star[(0,0)]**2,S_star+mx_star[(0,0)]**2])

    j +=1

mean = mx_star_mat[1:len(mx_star_mat)]

variance = S_star_mat[1:len(S_star_mat)]

mean1 = m_star_mat[1:len(m_star_mat)]

variance1 = Q_star_mat[1:len(Q_star_mat)]



xronos = np.arange(step,n+step,step)
#plt.ion()
plt.figure(5)
plt.plot(time1,trace)
##plt.plot(time1,ou,'r.')
plt.plot(xronos,obs,'ro')

plt.plot(t_all,y_all[:,0],'k')
plt.plot(t_all,y_all[:,0]-np.sqrt(y_all[:,1]),'g')
plt.plot(t_all,y_all[:,0]+np.sqrt(y_all[:,1]),'g')
##plt.plot(xronos,real_obs,'ro')
#plt.figure(2)
plt.title('KF2-OU')
#plt.plot(xronos,mean1)
#plt.xlim(0,2.5)
#plt.plot(xronos,mean-np.sqrt(variance),'g')
#plt.plot(xronos,mean+np.sqrt(variance),'g')
#plt.plot(xronos,mean)

#plt.ylim(-2,12)
plt.show()

#plt.figure(3)
#plt.plot(xronos,obs,'ro')
#plt.plot(xronos,mean1)
#plt.plot(xronos,mean1-np.sqrt(variance1),'g')
#plt.plot(xronos,mean1+np.sqrt(variance1),'g')
#plt.show()

##xronos = np.arange(2,26,2)
##plt.ion()
##plt.figure(3)
###plb.subplot(211)
##plt.plot(t_all,y_all[:,0])
###plb.pause(0.0001)  
###plt.plot(xronos,obs,'ro')
###plt.step(t_all,state[1:state.size+1,0])
##plt.xlim((0,24))
###plb.subplot(212)
##plt.figure(4)
##plt.plot(t_all,y_all[:,1])
###plb.pause(0.0001)  
###plt.step(t_all,state[1:state.size+1,1])
##plt.xlim((0,24))
###plb.show(block=True)
##plt.draw()
####plt.plot(t_all,y_all[:,0])
####plt.xlim((0,45))
####plt.show()
####plt.plot(t_all,y_all[:,1])
####plt.xlim((0,45))
####plt.show()
