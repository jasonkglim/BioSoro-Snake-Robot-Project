import numpy as np
import math
import matplotlib.pyplot as plt

# robot parameters
mu_t = 0.01
mu_n = 0.5
n = 3
link_mass = 1
link_length = 1
lG = 0.3 # distance to center of gravity of link
masses = np.repeat(link_mass, n)
lengths = np.repeat(link_length, n)
Glengths = np.repeat(lG, n)
L = np.sum(lengths)
mass_total = np.sum(masses)

# motion parameters
alpha_0 = np.pi / 3.
Kn = 1.
vel_s = 1.
phi_offset = 0
s = 0
init_pos = np.asmatrix(np.zeros(shape=(2,1)))
init_vel = np.asmatrix(np.zeros(shape=(2,1)))
joint_pos = np.array([np.asmatrix(np.zeros(shape=(2, 1))) for i in range(n+1)])
joint_vel = np.array([np.asmatrix(np.zeros(shape=(2, 1))) for i in range(n+1)])
joint_acl = np.array([np.asmatrix(np.zeros(shape=(2, 1))) for i in range(n+1)])
heading = 0
Bl = 2*Kn*np.pi/L
Bn = 2*Kn*np.pi/n

  
# initialize all matrices and vectors
# const matrices:
D = np.asmatrix(np.zeros(shape=(n, n-1)))
for i in range(n):
  for j in range(n-1):
    if i == j:
      D[i, j] = -1
    elif i == j+1:
      D[i, j] = 1
      
e = np.asmatrix(np.ones(shape=(n, 1)))

E = np.asmatrix(np.zeros(shape=(n, n-1)))
for i in range(n):
  for j in range(n-1):
    if i > j and i > 0:
      E[i, j] = 1
    else:
      E[i, j] = 0


M = np.asmatrix(np.ones(shape=(n, n)))

M_0 = np.asmatrix(np.zeros(shape=(n, 2)))

m_0 =  mass_total * np.asmatrix(np.diag(np.ones(2)))

m = np.asmatrix(np.zeros(shape=(2, n)))

f_f = np.asmatrix(np.ones(shape=(2, 1)))


f_0 = f_f

Phi = np.asmatrix(np.zeros(shape=(n, 1)))
Phidot = np.asmatrix(np.zeros(shape=(n, 1)))
Phiddot = np.asmatrix(np.zeros(shape=(n, 1)))

Phi[0] = alpha_0 * np.cos((2*Kn*np.pi/L) * s)  + phi_offset

Theta = np.asmatrix(np.zeros(shape=(n-1, 1)))
Thetadot = np.asmatrix(np.zeros(shape=(n-1, 1)))
Thetaddot = np.asmatrix(np.zeros(shape=(n-1, 1)))
for i in range(1, n):
  Theta[i-1] = -2*alpha_0*np.sin(Kn*np.pi/n) * np.sin(
    (2*Kn*np.pi/L) * s + (2*Kn*np.pi/n) * i ) + heading


# print "Theta: "
# print Theta

# print "Phi_0: ", Phi[0]

Phi = E * Theta + e * Phi[0]

# print "Phi: "
# print Phi

# print "a= ", np.sum([lengths[j]*np.cos(Phi[j]) for j in range(1)])
# print "b= ", np.sum([lengths[j]*np.sin(Phi[j]) for j in range(1)])

# print lengths[0]*np.cos(Phi[0])
# print lengths[0]*np.sin(Phi[0])

x_i = np.zeros(n+1)
y_i = np.zeros(n+1)

# Calculate joint positions
for i in range(n+1):
  a = np.sum([lengths[j]*np.cos(Phi[j]) for j in range(i)])
  b = np.sum([lengths[j]*np.sin(Phi[j]) for j in range(i)])
  joint_pos[i] = init_pos + np.asmatrix([[a], [b]])
  x_i[i] = joint_pos[i, 0]
  y_i[i] = joint_pos[i, 1]


print joint_pos
# print x_i
# print y_i


# plt.plot(x_i, y_i)
# plt.ylabel('y, meters')
# plt.xlabel('x, mmeters')
# plt.axis('equal')
# plt.grid()
# plt.show()



# torque inertial terms, mbar, M_0 and M
mbar = np.zeros(n)
for i in range(n):
  mbar[i] = masses[i] * Glengths[i] + lengths[i] * sum(masses[(i+1): n])

for i in range(n):
  M_0[i, 0] = -1*mbar[i] * np.sin(Phi[i])
  M_0[i, 1] = mbar[i] * np.cos(Phi[i])



x = init_pos.reshape(2)
y = np.array([1, 2])
# y = np.array([joint_pos[0, 0], joint_pos[0, 1]]).reshape(2)
print x
print y
print np.inner(x, y), np.sign(np.inner(x, y))




period = (L / vel_s)
tstep = period / 100
# Simulation starts here

# for t in range(0, period * 2, tstep):

  # s(t)
  s = vel_s * t
  
  # Calculate initial vectors
  Phi[0] = alpha_0 * np.cos((2*Kn*np.pi/L) * s)  + phi_offset  

  for i in range(1, n):
    Theta[i-1] = -2*alpha_0*np.sin(Bn/2) * np.sin(Bl* s + Bn*i) + heading

  Phi = E * Theta + e * Phi[0]

  Phidot[0] = -1*alpha_0*Bl* vel_s * np.sin(Bl * vel_s * t)
  init_vel = vel_s * np.asmatrix([[np.cos(Phi[0])], [np.sin(Phi[0])]])

  for i in range(1, n):
    Thetadot[i-1] = (-4*alpha_0*Kn*np.pi/L)*np.sin(Kn*np.pi/n) * np.cos(Bl*s + Bn*i)

  Phidot = E*Thetadot + e*Phidot[0]

  for i in range(n+1):
    a = np.sum([lengths[j]*np.cos(Phi[j]) for j in range(i)])
    b = np.sum([lengths[j]*np.sin(Phi[j]) for j in range(i)])
    joint_pos[i] = init_pos + np.asmatrix([[a], [b]])

    a = np.sum([-1*lengths[j]*np.sin(Phi[j])*Phidot[j] for j in range(i)])
    b = np.sum([lengths[j]*np.cos(Phi[j])*Phidot[j] for j in range(i)])
    joint_vel[i] = init_vel + np.asmatrix([[a], [b]])

    f_posdot[i] = joint_vel + int(flengths[i]*Phidot[i]) * np.asmatrix([[np.cos(Phi[i])], [np.sin(Phi[i])]])
    u_t[i] = np.array([np.cos(Phi[i]), np.sin(Phi[i])])
    u_n[i] = np.array([-1*np.sin(Phi[i]), np.cos(Phi[i])])
    sign_t[i] = int(np.sign(np.inner(f_posdot[i].reshape(2), u_t[i])))
    sign_n[i] = int(np.sign(np.inner(f_posdot[i].reshape(2), u_n[i])))

    
