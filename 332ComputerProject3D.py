import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

#choose the mode of the standing waves for each direction
#I chose them such that the frequency is nearest to 2.4GHz which is what a regular microwave operates at
n_x=3
n_y=1
n_z=3

#choose x length of microwave
L_x=.35
L_y=.33 #L_x*0.9
L_z=.21 #L_x*0.7

k_x=(n_x*np.pi)/L_x
k_y=(n_y*np.pi)/L_y
k_z=(n_z*np.pi)/L_z

#mu pretty much cancels itself out with the speed of light
mu = 4*np.pi*10**(-7)
c=3*10**8
w = c*np.sqrt(k_x**2 + k_y**2 + k_z**2)


x = np.linspace(0,L_x,1000)
y = np.linspace(0,L_y,1000)
z = np.linspace(0,L_z,1000)
xy, yx = np.meshgrid(x,y,sparse=True)
xz, zx = np.meshgrid(x,z,sparse=True)
yz, zy = np.meshgrid(y,z,sparse=True)


def s_z(t,z,x,y):
    s_z = (1/mu)*(1/w)*np.cos(w*t)*np.sin(w*t)*np.sin(k_z*z)*np.cos(k_z*z)*((np.cos(k_x*x)**2)*(np.cos(k_y*y)**2)*(k_x+k_y-(2*k_z)) + ((np.cos(k_x*x)**2)*(k_z-k_x)) + ((np.cos(k_y*y)**2)*(k_z-k_y)))
    return s_z

def s_y(t,z,x,y):
    s_y = (1/mu)*(1/w)*np.cos(w*t)*np.sin(w*t)*np.sin(k_y*y)*np.cos(k_y*y)*((np.cos(k_x*x)**2)*(np.cos(k_z*z)**2)*(k_x+k_z-(2*k_y)) + ((np.cos(k_z*z)**2)*(k_y-k_z)) + ((np.cos(k_x*x)**2)*(k_y-k_x)))
    return s_y

def s_x(t,z,x,y):
    s_x = (1/mu)*(1/w)*np.cos(w*t)*np.sin(w*t)*np.sin(k_x*x)*np.cos(k_x*x)*((np.cos(k_y*y)**2)*(np.cos(k_z*z)**2)*(k_y+k_z-(2*k_x)) + ((np.cos(k_y*y)**2)*(k_x-k_y)) + ((np.cos(k_z*z)**2)*(k_x-k_z)))
    return s_x

def poyntingMag(t):
    S=np.sqrt(s_z(t,z,xy,yx)**2 + s_y(t,z,xy,yx)**2 + s_x(t,z,xy,yx)**2)
    return S

def Surface_Area_Calc_z(S,xvals,yvals,x_i, x_f, y_i, y_f, z):
    power = S(np.pi/8,z,xvals,yvals)[x_i:x_f][y_i:y_f]
    return power

def Surface_Area_Calc_y(S,xvals,zvals,x_i,x_f,z_i,z_f,y):
    power = S(np.pi/8,zvals,xvals,y)[x_i:x_f][z_i:z_f]
    return power

def Surface_Area_Calc_x(S,yvals,zvals,y_i,y_f,z_i,z_f,x):
    power = S(np.pi/8,zvals,x,yvals)[y_i:y_f][z_i:z_f]
    return power
#I know the length of the lists are all 50 so i can choose the spots to include in the boundary of the power calc
#I know the avg power happens at t=pi/8 for the time dependence sin(t)cos(t)


#Will assume the object in the food in the microwave is the size of a hot pocket which I am estimating as the size of my phone
#I can not find the size anywhere (phone dimensions -> 0.1589 X 0.0736 X 0.0084 m) {don't microwave your phone}


def power_to_object(x1,x2,y1,y2,z1,z2):
    power = 0
    #This is pretty much a bunch of dot products
    #Taking x component of poynting vector and multiplying it by surfaces with direction x-hat etc...

    #I dont include a bottom surface because it is usually on a plate and cant receive any reflected waves if I dont consider transmission
    top_surf = (Surface_Area_Calc_z(s_z,xy,yx,x1,x2,y1,y2,z2))
    power += np.sum(abs(top_surf))
    back_surf = (Surface_Area_Calc_y(s_y,xz,zx,x1,x2,z1,z2,y1))
    power += np.sum(abs(back_surf))
    front_surf = (Surface_Area_Calc_y(s_y,xz,zx,x1,x2,z1,z2,y2))
    power += np.sum(abs(front_surf))
    #The left and right side add very little to power because they are the smallest sides
    left_surf = (Surface_Area_Calc_x(s_x,yz,zy,y1,y2,z1,z2,x1))
    power += np.sum(abs(left_surf))
    right_surf = (Surface_Area_Calc_x(s_x,yz,zy,y1,y2,z2,z2,x2))
    power += np.sum(abs(right_surf))

    return power


#x1=.09555->273, x2=.25445->727, y1=.1282->384.5, y2=.2018->611.5, z1=0->0, z2=0.0084->40
#can run this function to put the hot pocket at certain coordinates
print('Average power delivered to a centered hotpocket:',round(power_to_object(273,727,385,812,0,40),4),'W')#this describes a centered hotpocket with the longest side along the widest microwave dimension

#I decided not to use thee because there are too many things run
#with the hastag removed different positions specifically can be shown
#power_to_object(23,477,385,812,0,40)#this is halfway away from the center and left
#power_to_object(0,454,385,812,0,40)#this describes the hot pocket being on the 'left' end of the microwave 
#power_to_object(523,977,385,812,0,40)#halfway between the middle and right side
#power_to_object(546,1000,385,812,0,40)#this is at the 'right' side of the microwave

power_list_x = []
xlist = []
x_pos = np.arange(0,546,10)
for x in x_pos:
    xlist.append(x+227)
    power_list_x.append(power_to_object(0+x,454+x,385,812,0,40))

plt.plot(xlist,power_list_x)
plt.xlabel('position of center of hotpocket changing only in x direction')
plt.ylabel('Average Power [$W$]')
plt.show()

power_list_y = []
ylist = []
y_pos = np.arange(0,573,10)
for y in y_pos:
    ylist.append(y+213)
    power_list_y.append(power_to_object(273,727,0+y,427+y,0,40))

plt.plot(ylist,power_list_y)
plt.xlabel('position of center of hotpocket changing only in y direction')
plt.ylabel('Average Power [$W$]')
plt.show()

power_list_z = []
zlist = []
z_pos = np.arange(0,960,10)
for z in z_pos:
    zlist.append(z+20)
    power_list_z.append(power_to_object(273,727,385,812,0+z,40+z))

plt.plot(zlist,power_list_z)
plt.xlabel('position of center of hotpocket changing only in z direction')
plt.ylabel('Average Power [$W$]')
plt.show()

indexx = power_list_x.index(max(power_list_x))
indexy = power_list_y.index(max(power_list_y))
indexz = power_list_z.index(max(power_list_z))

print('The best position to put the hotpocket is at', round((.35*indexx/1000),5),'cm from the left end', round((.33*indexy/1000),5),'cm from the front', round((.21*indexz/1000),5), 'cm from the bottom')

print('Maximum average power delivered to a hotpocket:',round(power_to_object(indexx, indexx+454, indexy, indexy+426, indexz, indexz+40),4),'W')
#I don't think this is the true maximum but looking at the max for each direction as it passes through the center this is true
#I dont think i have the power to compute the maximum possible from the entire microwave
#Mostly think its important to show that at this 'fake' maximum it is better than the centered hotpocket