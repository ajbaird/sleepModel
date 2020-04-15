import matplotlib.pyplot as plt
import numpy as np 
from scipy.integrate import solve_ivp

# system of equations 
def circadianProcess (time):
    params = [5.0, 2.5, 2.0]
    mc = params[0]
    A = params[1]
    sleepTime = params[2]

    #scaling this equation to get it in units of hours
    c = mc - (0.8*A)*np.sin((2*np.pi/24.0)*time*(0.5*(sleepTime/24.0)))

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
def eulerStep(BD, rbt, rwt, pw, pb1, dt, time, index):

    #hard code these for now:


    ct = circadianProcess(time)
    xb = ct * (BD[index] / (1 + BD[index]**2))
    
    #update the next time step
    BD[index+1] = BD[index] + dt*(pw*rwt + pb1*rbt*BD[index] - rbt*xb)


def computeDebt(t, dt, sleepAmount, paramrbt, paramrwt, parampw, parampb1):

    temp = 0
    BD = np.zeros(len(t) + 1)
    #set initial value
    BD[0] = 0.0
    sleepstart = 24.0 - sleepAmount

    for i in t: 
        if sleepstart < i < 24.0 or (48 - sleepAmount) < i < 48.0 or (72 - sleepAmount) < i < 72.0:
         #scale rbt and rwt during sleep periods:
            rbt = 0.28 * paramrbt
            rwt = 0.06 * paramrwt
            pw = 0.13 * parampw
            pb1 = 1.7 * parampb1
        else: 
            rbt = 0.28
            rwt = 0.06
            pw = 0.13 
            pb1 = 1.7

        #take an euler step
        eulerStep(BD, rbt, rwt, pw, pb1, dt, i, temp)
        temp += 1 

    #return solution
    return BD

#main: 
if __name__ == "__main__":
    
#grab a domain (24 hour window here)
    t, dt = np.linspace(0, 72, 10000, retstep = True) 

    sleepamount = 2.0
    paramrbt = 1.1#1.00
    paramrwt = 0.01#1.0
    parampw = 1.0#1.0
    parampb1 = 1.0
    sleepstart = 24.0 - sleepamount
    bioDebtNormal = computeDebt(t, dt, sleepamount, paramrbt, paramrwt, parampw, parampb1)
    bioDebtSleepy = computeDebt(t, dt, sleepamount/2.0, paramrbt, paramrwt, parampw, parampb1)
    bioDebtAwake =  computeDebt(t, dt, sleepamount * 2.0, paramrbt, paramrwt, parampw, parampb1)

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

    fig, ax = plt.subplots(3, sharex=True, sharey=True)

    ax[0].plot(t, bioDebtNormal[:-1])
    ax[0].set_title('Debt')

    ax[1].plot(t, bioDebtSleepy[:-1])
    ax[1].set_title('Debt sleep half')

    ax[2].plot(t, bioDebtAwake[:-1])
    ax[2].set_title('Debt sleep double')

    plt.tight_layout()
    #plt.savefig('HighA2HoursRateLimitRW0p1.png', dpi=300 )

    #plt.savefig('lowA8HoursRateLimitRW0p1.png', dpi=300 )
    #plt.savefig('highA2HoursRateLimitRW0p01.png', dpi=300 )


    plt.show()