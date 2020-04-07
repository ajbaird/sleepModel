import matplotlib.pyplot as plt
import numpy as np 
from scipy.integrate import solve_ivp

# system of equations 
def circadianProcess (time):
    params = [5.0, 2.5, 8.0]
    mc = params[0]
    A = params[1]
    sleepTime = params[2]

    #scaling this equation to get it in units of hours
    c = mc - A*np.sin((2*np.pi/24.0)*(time - 0.25 + 0.5*(sleepTime/24.0)))

    return c

def sleepModel(t, x, params):

    
    I = x[0]
    R = x[1]
    D = x[2]

    ct = circadianProcess(t)
    #print(ct)
    xb = ct * (D / (1 + D**2))
    pb = params[3] * D
    D = R - I
    dIdt = xb*params[2]
    dRdt = params[0]*params[1] + pb * params[2]
    dDdt = params[0]*params[1] + (pb - xb)*params[2]

    return [dIdt, dRdt, dDdt]

#custom evaluations in order to change parameters to reflect metabolic rate
def eulerStep(BD, rbt, rwt, dt, time, index):

    #hard code these for now:
    pw = 0.13
    pb1 = 1.7

    ct = circadianProcess(time)
    xb = ct * (BD[index] / (1 + BD[index]**2))
    
    #update the next time step
    BD[index+1] = BD[index] + dt*(pw*rwt + pb1*rbt*BD[index] - rbt*xb)


def computeDebt(t, dt, sleepAmount):

    temp = 0
    BD = np.zeros(len(t) + 1)
    #set initial value
    BD[0] = 0.0
    sleepstart = 24.0 - sleepAmount

    for i in t: 
        if sleepstart < i < 24.0 or sleepstart*2.0 < i < 48.0 or sleepstart*3.0 < i < 48.0:
         #scale rbt and rwt during sleep periods:
            rwt = 0.06 * 1.0
            rbt = 0.28 * 0.0010
        else: 
            rbt = 0.28
            rwt = 0.06

    #take an euler step
        eulerStep(BD, rbt, rwt, dt, i, temp)
        temp += 1 

    #return solution
    return BD

#main: 
if __name__ == "__main__":
    
#grab a domain (24 hour window here)
    t, dt = np.linspace(0, 72, 10000, retstep = True) 

    sleepamount = 4.0
    sleepstart = 24.0 - sleepamount
    bioDebt = computeDebt(t, dt, sleepamount)

    #params stored in a vector, passed to model, pw, rw, rb, pb1: 
    params = [0.13, 0.06, 0.28, 1.7]

    #intial conditions I, R, D:  
    x0 = [0.1, 0.2, 0.0]

    #compute bd via ivp solver (for testing)
    sol = solve_ivp(lambda t, x: sleepModel(t, x, params), [t[0], t[-1]], x0,  t_eval = t )


#collect terms 
I = sol.y[0]
R = sol.y[1]
D = sol.y[2]
c = circadianProcess(t)

fig, ax = plt.subplots(2,2)

ax[0,0].plot(t, I)
ax[0,0].set_title('Investments')

ax[0,1].plot(t, R)
ax[0,1].set_title('Requirements')

ax[1,0].plot(t, D)
ax[1,0].set_title('Debt')

#ax[1,1].plot(t, c)
#ax[1,1].set_title('circadian process')

ax[1,1].plot(t, BD[:-1])
ax[1,1].set_title('Debt New')

plt.show()