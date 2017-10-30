'''
The first part of the program sets the initial parameters for an object thrown at some angle above the horizont, including the air resistance dragging effect.
Then, it solves the differential equation m(dV/dT)=mg-cV, describing the motion ( the numerical solution).
The second part computes the analytical solution and plots both of them in a x-y plot. 
The third part creates an animation of the motion and shows the landing point. 

The new things used are odeint, plots with title, grid and labeled arrow, and animation.
'''
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.animation as animation


c=0.7 #Sets the dragging term.
m=0.1 #Sets object's mass.
g=9.81 #Sets gravitational acceleration.
Vx0=10 #Sets initial speed along the x-axis.
Vy0=10 #Sets initial speed along the y-axis.
thetha=math.pi/4 #Sets the angle of throwing.
V0=np.sqrt(Vx0**2+Vy0**2) #Calculates the initial velocity.

t0 = 0.0  #  Sets the initial time
tmax=2.0  #  Sets the final final time
steps=100  # Sets the number of time step
tLF = np.linspace(t0, tmax, steps+1)  # Creates a 1-D array of time values

y0LF = [0.0, Vx0, 0.0, Vy0]  # Creates an array with the initial condition for x-position, velocity along x, y-position, velocity along y.

def deriv(yLF,tF): #Creates the function which computes the derivatives 
    Vx = yLF[1]   # Identifies the velocity along x axis
    Vy=yLF[3]   #Identifies the velocity along y axis
    return [Vx, -c*Vx/m, Vy, -g-c*Vy/m]   # The function outputs [dx/dt, dVx/dt, dy/dt, dVy/dt]

yM = odeint(deriv, y0LF, tLF)  #  The 4-D array containing the solution for the differential equation

plt.plot(yM[:,0],yM[:,2],'.',label='Numerical solution') #Plots y over x numerically.

#Analytical Solution:
VT=m*g/c #Calculates the terminal velocity
anal_x=((V0*VT)/g)*np.cos(thetha)*(1-np.exp(-g*tLF/VT)) #calculates dx/dt using the analytical solution
anal_y=(VT/g)*(V0*np.sin(thetha)+VT)*(1-np.exp(-g*tLF/VT))-VT*tLF #calculates dy/dt using the analytical solution

plt.plot(anal_x,anal_y,label='Analytical solution') #Plots y over x analytically.
plt.grid()
plt.xlabel('Horizontal axis [meters]')
plt.ylabel('Vertical axis [meters]')
plt.legend()
plt.title('Projectile motion with air resistance.')
plt.savefig('ballistic.pdf')

#Computing answer of the questions:

i=np.abs(yM[1:steps,2]).argmin() #calculates the point from the y-array closest to 0 (after the initial point)and assigns it to i (=impact).

D=yM[i+1,0] #Computes the distance to the point of impact.

H=np.amax(yM[:,2]) #Finds the max value of the array anal_y, which represents the highest point of the trajectory and assignes it to H.

TF=tLF[i+1]

Vxi=yM[i+1,1]
Vyi=yM[i+1,3]
Vi=np.sqrt(Vxi**2+Vyi**2)

#Creates the file .txt and saves the answers of the questions
f=open('ballistic.txt','w')
f.write(' The distance to the point of impact is {}.\n The highest point of the trajectory is at {} meters above ground.\n The time for flight is {} s. \n The impact velocity is {} m/s.'.format(str(D),str(H),str(TF),str(Vi)))
f.close()


#Extra: animation of the projectile motion.
fig, ax = plt.subplots() #sets the initial figure.
line, = ax.plot(anal_x,anal_y)   # Sets the values of 'line'

def animate(k): #Creates the animation function, which updates the data in each frame. 
    line.set_data(anal_x[:k],anal_y[:k])
    return line,

plt.axis([0.0, 2.5, 0.0, 2.5]) #Sets the range of the axis. 

ani = animation.FuncAnimation(fig, animate,100,interval=50,blit=True) #Creates the animation. 
plt.title('Animation of projectile motion with air resistance.')
#Adds an arrow pointing at the landing spot:
plt.annotate('landing point', xy=(D, 0), xytext=(1.8, 0.8), arrowprops=dict(facecolor='blue', shrink=0.02),)
plt.show()
