import numpy as np
import matplotlib.pyplot as plt

class Orbit():
  #                    int               array         float      array                   float
  def __init__(self, number_of_moons, inital_position, mass, velocity_at_inital_position, dt, eccentricity):
    self.moons = number_of_moons
    self.position = np.copy(inital_position)
    self.mass = mass
    self.velocity = np.copy(velocity_at_inital_position)
    self.dt = dt
    self.e = eccentricity

  def magnitude(self, r):
    return np.sqrt(r[0]**2 + r[1]**2)
  
  def finddr(self,r1,r2):
        return np.sqrt((r2[0]-r1[0])**2 + (r2[1]-r1[1])**2)
  
  def derivs(self, vars, other_object):
    r = vars[:2]
    v = vars[2:]
    
    mass_sun = 2.0E30 #kg
    
    xDeriv = v[0]
    yDeriv = v[1]

    vxDeriv = (-4 * (np.pi**2) *r[0] /(abs(self.magnitude(r))**3)) + (4 * (np.pi**2) * ((other_object.position[0] - r[0])* other_object.mass) / (mass_sun * abs(self.finddr(r,other_object.position))**3 ))
    vyDeriv = (-4 * (np.pi**2) *r[1] /(abs(self.magnitude(r))**3)) + (4 * (np.pi**2) * ((other_object.position[1] - r[1])* other_object.mass) / (mass_sun * abs(self.finddr(r,other_object.position))**3 ))

    return np.array([xDeriv,yDeriv,vxDeriv,vyDeriv])

  def orbit_path(self, other_object):
    vars = np.array([self.position[0],self.position[1],self.velocity[0],self.velocity[1]])
    
    k1 = self.dt * self.derivs(vars, other_object)
    k2 = self.dt * self.derivs(vars + 1/2 * k1, other_object)
    k3 = self.dt * self.derivs(vars + 1/2 * k2, other_object)
    k4 = self.dt * self.derivs(vars + k3, other_object)
    vars += 1/6 * (k1 + 2 * k2 + 2 * k3 + k4)
    #vars += k2
    self.position = vars[:2]
    self.velocity = vars[2:]

    return self.position
    

ee = 0.017
em = 0.0549

Erad = 1
Mrad = 0.0025

a = 0.0
b = 1
h = 1 / 10000
tpoints = np.arange(a,b,h)


earth = Orbit(1, np.array([1 * (1+ee) , 0]) ,6.0e24 , np.array([0 , np.sqrt((4* np.pi**2)*(1-ee) / (Erad*(1+ee)))]) , h , 0.017)

moon = Orbit(0, np.array([earth.position[0] + (0.0025) * (1+em) , 0]) , 7.3e22 , np.array([0 , np.sqrt((4* np.pi**2 * (earth.mass/2.0E30) * (1-em))/(Mrad*(1+em))) + earth.velocity[1]]) , h , 0.0549)


earth_pathxm =[]
earth_pathym = []
moon_pathxm = []
moon_pathym = []

for t in tpoints: # 1 year
  em = earth.orbit_path(moon)
  mm = moon.orbit_path(earth)

  earth_pathxm.append(em[0])
  earth_pathym.append(em[1])

  moon_pathxm.append(mm[0])
  moon_pathym.append(mm[1])



fig,(ax1) = plt.subplots(1,1,figsize=(10,10))

ax1.scatter(0,0, color = 'orange', linewidth = 15, label ="Sun")
ax1.scatter(earth_pathxm,earth_pathym,label='Earth', linewidths=10)
ax1.scatter(moon_pathxm,moon_pathym,label='Moon', linewidths=9, s = .1)
plt.xlabel("X position (AU)")
plt.ylabel("Y position (AU)")
ax1.legend()
plt.show()

fig,(ax1) = plt.subplots(1,1,figsize=(10,10))

ax1.scatter(0,0, color = 'orange', linewidth = 15, label ="Sun")
ax1.scatter(earth_pathxm,earth_pathym,label='Earth', linewidths=10)
ax1.scatter(moon_pathxm,moon_pathym,label='Moon', linewidths=9, s = .1)
plt.xlabel("X position (AU)")
plt.ylabel("Y position (AU)")
plt.xlim(0,1)
plt.ylim(0,1)
ax1.legend()
plt.show()