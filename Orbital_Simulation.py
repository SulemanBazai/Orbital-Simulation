"""
Goal/Purpose:  Determining the orbit circumference  and perigee precession
               using adaptive midpoint integration and fine scanning
               with post Newtonian correction to force of gravity

Author: Suleman Bazai

Date: October 18, 2018
    
Physics 298 owl, University of Illinois at Urbana-Champaign, 2018
"""

# import libraries
import numpy as np
import matplotlib.pyplot as plt 
import time


# Define and initialize variables
start_time = time.time()

 # Newton's constant
G = 6.67408E-11

earth_radius = 6.371E6 #in meters

#earth_mass in kg
m_earth = 5.972E24

#post newtonian factor
post_newtonian_factor = 1000.0

#function to find x and y accelerations given the current position (x,y) and the mass of the earth
def force_function(x,y):
    ax = (-G*m_earth*x)/((x**2+y**2)**(3/2))\
        -(post_newtonian_factor*x*3*8.87e-3*(6.471e6*10800)**2)\
        /(2*(x**2+y**2)**(5/2))
    ay = (-G*m_earth*y)/((x**2+y**2)**(3/2)) \
        -(post_newtonian_factor*y*3*8.87e-3*(6.471e6*10800)**2)\
        /(2*(x**2+y**2)**(5/2))
    return ax,ay
    
#initial conditions
input_v0 = 1.08e4
t = 0
x = 6.471e6
y = 0
vy = input_v0
vx = 0
r = 6.471e6
dt = .1
dx = 0
dy = 0
i=0
dt_counter = 0
bins = 0
circumference = 0
circumference_counter = 0

#to generate big enough arrays for time
tmax = 3600000

#initialize position and velocity arrays
x_array = np.empty(tmax)
y_array = np.empty(tmax)
vx_array = np.empty(tmax)
vy_array = np.empty(tmax)
t_array = np.empty(tmax)
r_array = np.zeros(tmax)
dx_array = np.empty(tmax)
dy_array = np.empty(tmax)
dt_array = np.empty(tmax)
sigma_array = np.zeros(tmax)
d_circumference_array = np.zeros(tmax)

#max time we are interested in in seconds
tmax = 2*24*3600



#############################################################################
#loop to iterate through the motion and position of the satellite in its orbit
#using the adaptive midpoint integration techqniue with fine scanning
# in order to find the precession during an orbit and the orbit circumference
#############################################################################
while t<tmax:
    #update arrays at the start of the time step
    x_array[i] = x
    y_array[i] = y
    vx_array[i] = vx
    vy_array[i] = vy
    t_array[i] = t
    dx_array[i] = dx
    dy_array[i] = dy
    dt_array[i] = dt
    r_array[i] = r
    
    # calculate accelerations at start of this interval
    ax,ay = force_function(x ,y)
    
    # calculate x_midpoint and v_midpoint (i.e at time step of dt/2)
    xmid = x + vx * dt / 2
    ymid = y + vy * dt / 2
    vxmid = vx + ax * dt / 2
    vymid = vy + ay * dt / 2
    
    # calculate accelerations at midpoint of this interval
    axmid,aymid = force_function(xmid,ymid)
    
    # calculate approximate changes to x,y, vx, vy over the full interval using the midpooint velocities and accerlerations
    dx = vxmid*dt
    dy = vymid*dt
    dvx = axmid*dt
    dvy = aymid*dt
    
    # update position, velocity, time, index, and change dt by appropriate amount
    x = x + dx
    y = y + dy
    r = np.sqrt(x**2+y**2)
    vx = vx + dvx
    vy = vy + dvy
    t = t + dt 
    dt_counter+=1
    d_circumference_array[i] = np.sqrt(dx**2+dy**2)
    
    #conditoinal statements to determine whether or not to change the time step dt for the next step of the integration
    if -25000<y<13000 and x>0 and i>1000:
        dt = 1e-5
        bins+=1
    elif dt_counter>2:
        sigma_array[i] = (1/81)*np.sqrt(((2*dx_array[i]-dx_array[i-1]-dx)**2)\
                   +((2*dy_array[i]-dy_array[i-1]-dy)**2))
        if sigma_array[i]>1e-6:
            dt = dt/2
            dt_counter = 0
        elif sigma_array[i]<1e-9:
            dt = dt*2
            dt_counter = 0
        else:
            dt = dt
    else:
        sigma_array[i]=sigma_array[i-1]
              
    i += 1


indexmin = np.argmin(r_array[13000:i])+13000
#loop through change in circumference through each iteration to find orbit 
#circumference
while circumference_counter<indexmin:
    circumference += d_circumference_array[circumference_counter]
    circumference_counter +=1
    
rmax = np.max(r_array[:i])
rmin = np.min(r_array[:i])
closure_miss_distance = np.sqrt((x_array[indexmin]-x_array[0])**2\
                                +(y_array[indexmin]-y_array[0])**2)
elapsed_time = time.time()-start_time

print("Post-Newtonian correction scale factor =",post_newtonian_factor)
print("Maximum orbital radius (km):",rmax/1000)
print("Minimum orbital radius (km):",rmin/1000)
print("Number of bins in the fine-grained scan:",bins)
print("Index_perigee:",indexmin, "    value:",r_array[indexmin])
print("Closure miss distance:", closure_miss_distance)
print("Orbit circumference (km):", circumference/1000)
print("Apparent perigee precession (microradians), one orbit:"\
      ,(y_array[indexmin])/(x_array[indexmin])/1e-6)
print("Elapsed running time:",elapsed_time,"seconds")




