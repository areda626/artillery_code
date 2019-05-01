# start from Example 8.5 program odesim.py
import numpy as np
import matplotlib.pyplot as plt

# constants
g = 9.81 # m s^-2
g_mars  = 3.71	# m s^-2
m     = 1.0	# kg
rho_mars   = 0.20	# kg m^-3
rho = 1.22# kg m^-3
C     = 0.47	# unitless
R     = 0.08    # m
h     = 0.001   # seconds
theta = 30.0*(np.pi/180) # radians
v0    = 100.0	# m s^-1
const_mars = (rho_mars*C*np.pi*R**2)/(2.0*m)
const = (rho*C*np.pi*R**2)/(2.0*m) 

# define the equations of motion
def f_mars(r,const):
    x   = r[0]
    y   = r[1]
    vx  = r[2]
    vy  = r[3]
    fx  = vx
    fy  = vy
    fvx = -const*vx*np.sqrt(vx**2+vy**2)
    fvy = -g_mars-const*vy*np.sqrt(vx**2+vy**2)
    return np.array([fx,fy,fvx,fvy],float)

# define the equations of motion
def f(r,const):
    x   = r[0]
    y   = r[1]
    vx  = r[2]
    vy  = r[3]
    fx  = vx
    fy  = vy
    fvx = -const*vx*np.sqrt(vx**2+vy**2)
    fvy = -g-const*vy*np.sqrt(vx**2+vy**2)
    return np.array([fx,fy,fvx,fvy],float)
# containers for output
r = np.array([0.0,0.0,v0*np.cos(theta),v0*np.sin(theta)],float)
xpoints = []
ypoints = []

# use fourth-order Runge-Kutta
while r[1]>=0:
    k1 = h*f(r,const)
    k2 = h*f(r+0.5*k1,const)
    k3 = h*f(r+0.5*k2,const)
    k4 = h*f(r+k3,const)
    r += (k1+2*k2+2*k3+k4)/6
    xpoints.append(r[0])
    ypoints.append(r[1])

# make plot for part (b)
p1 = plt.figure(1)
plt.plot(xpoints,ypoints)
plt.xlabel("x [m]")
plt.ylabel('y [m]')
p1.show()
p1.savefig('earth.png')

r = np.array([0.0,0.0,v0*np.cos(theta),v0*np.sin(theta)],float)
xpoints = []
ypoints = []
# use fourth-order Runge-Kutta
while r[1]>=0:
    k1 = h*f_mars(r,const_mars)
    k2 = h*f_mars(r+0.5*k1,const_mars)
    k3 = h*f_mars(r+0.5*k2,const_mars)
    k4 = h*f_mars(r+k3,const_mars)
    r += (k1+2*k2+2*k3+k4)/6
    xpoints.append(r[0])
    ypoints.append(r[1])

# make part b plot for mars
p2 = plt.figure(2)
plt.plot(xpoints,ypoints)
plt.xlabel("x [m]")
plt.ylabel('y [m]')
p2.show()
p2.savefig('mars.png')

# try different values of m
p3 = plt.figure(3)
for m in [0.25,0.5,1,2,4]:
    const = (rho*C*np.pi*R**2)/(2.0*m)
    r = np.array([0.0,0.0,v0*np.cos(theta),v0*np.sin(theta)],float)
    xpoints = []
    ypoints = []

    # use fourth-order Runge-Kutta
    while r[1]>=0:
        k1 = h*f(r,const)
        k2 = h*f(r+0.5*k1,const)
        k3 = h*f(r+0.5*k2,const)
        k4 = h*f(r+k3,const)
        r += (k1+2*k2+2*k3+k4)/6
        xpoints.append(r[0])
        ypoints.append(r[1])

    plt.plot(xpoints,ypoints,label='m = '+str(m)+' kg')

plt.xlabel("x [m]")
plt.ylabel('y [m]')
plt.legend()
p3.show()
p3.savefig('earth_multi.png')
# try different values of m for mars
p4= plt.figure(4)
for m in [0.25,0.5,1,2,4]:
    const = (rho_mars*C*np.pi*R**2)/(2.0*m)
    r = np.array([0.0,0.0,v0*np.cos(theta),v0*np.sin(theta)],float)
    xpoints = []
    ypoints = []

    # use fourth-order Runge-Kutta
    while r[1]>=0:
        k1 = h*f_mars(r,const)
        k2 = h*f_mars(r+0.5*k1,const)
        k3 = h*f_mars(r+0.5*k2,const)
        k4 = h*f_mars(r+k3,const)
        r += (k1+2*k2+2*k3+k4)/6
        xpoints.append(r[0])
        ypoints.append(r[1])

    plt.plot(xpoints,ypoints,label='m = '+str(m)+' kg')

plt.xlabel("x [m]")
plt.ylabel('y [m]')
plt.legend()
p4.show()
p4.savefig('mars_multi.png')
