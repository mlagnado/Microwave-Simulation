import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

#choose the mode of thestanding waves for each direction
#these selected modes get the frequency closest to that which is expected of a microwave for the chosen size
n_x=3
n_y=1
n_z=3

#choosen dimensions of microwave 35x33x21 cm
L_x=.35
L_y=.33#L_x*0.9
L_z=.21#L_x*0.7

k_x=(n_x*np.pi)/L_x
k_y=(n_y*np.pi)/L_y
k_z=(n_z*np.pi)/L_z

#mu pretty much cancels itself out with the speed of light so set both of them equal to 1
mu = 4*np.pi*10**(-7)
c=3*10**8
w = c*np.sqrt(k_x**2 + k_y**2 + k_z**2)


x = np.linspace(0,L_x,1000)
y = np.linspace(0,L_y,1000)
xx, yy = np.meshgrid(x,y,sparse=True)


#Looking at the base of the microwave -> z=0
z=0


def s_z(t,z,xx,yy):
    s_z = (1/mu)*(1/w)*np.cos(w*t)*np.sin(w*t)*np.sin(k_z*z)*np.cos(k_z*z)*((np.cos(k_x*xx)**2)*(np.cos(k_y*yy)**2)*(k_x+k_y-(2*k_z)) + ((np.cos(k_x*xx)**2)*(k_z-k_x)) + ((np.cos(k_y*yy)**2)*(k_z-k_y)))
    return s_z

def s_y(t,z,xx,yy):
    s_y = (1/mu)*(1/w)*np.cos(w*t)*np.sin(w*t)*np.sin(k_y*yy)*np.cos(k_y*yy)*((np.cos(k_x*xx)**2)*(np.cos(k_z*z)**2)*(k_x+k_z-(2*k_y)) + ((np.cos(k_z*z)**2)*(k_y-k_z)) + ((np.cos(k_x*xx)**2)*(k_y-k_x)))
    return s_y

def s_x(t,z,xx,yy):
    s_x = (1/mu)*(1/w)*np.cos(w*t)*np.sin(w*t)*np.sin(k_x*xx)*np.cos(k_x*xx)*((np.cos(k_y*yy)**2)*(np.cos(k_z*z)**2)*(k_y+k_z-(2*k_x)) + ((np.cos(k_y*yy)**2)*(k_x-k_y)) + ((np.cos(k_z*z)**2)*(k_x-k_z)))
    return s_x

def poyntingMag(t):
    S=np.sqrt(s_z(t,z,xx,yy)**2 + s_y(t,z,xx,yy)**2 + s_x(t,z,xx,yy)**2)
    return S

def update(t):
    ax.clear()
    data = poyntingMag(t)
    ax.contourf(x,y,data, cmap='RdYlBu_r', levels = 10)#, vmin=min, vmax=max)
    ax.set_ylabel('Distance into the microwave [m]')
    ax.set_xlabel('Distance along the microwave [m]')
    ax.set_title(f"Power delivered to the floor of a microwave at Time: {t:.1f} seconds")
    return ax


time_steps = np.arange(0,(2*np.pi),(np.pi/30))

fig, ax = plt.subplots()
contour = ax.contourf(x,y,poyntingMag(np.pi/4), cmap='RdYlBu_r')
colorbar = fig.colorbar(contour)
colorbar.set_label('Power [$W$]')
ani = FuncAnimation(fig,update,frames=time_steps, interval=500, blit=False)
plt.show()
print('This is the frequency of the standing wave in the microwave in GHz:', w/(2*np.pi*10**9))