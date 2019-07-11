import numpy as np
import math
import matplotlib.pyplot as plt

# robot parameters
mu_t = 0.01
mu_n = 0.65
n = 8
link_mass = 1.
link_length = 1.
link_radius = 0.1
lG = 0.5 # distance to center of gravity of link
masses = np.repeat(link_mass, n)
lengths = np.repeat(link_length, n)
Glengths = np.repeat(lG, n)
L = np.sum(lengths)
mass_total = np.sum(masses)

# motion parameters
incline = 0 #(np.pi/180) * 10
alpha_0 = np.pi / 6.
Kn = 1
vel_s = 1.2
phi_offset = 0
heading = 0.1
Bl = 2*Kn*np.pi/L
Bn = 2*Kn*np.pi/n
period = (L / vel_s)
N_t = 100
tstep = period / N_t
g = np.asmatrix([[9.8*np.sin(incline)], [0]])

  
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


# position, velocity and acceleration lists of vectors
init_pos = np.asmatrix(np.zeros(shape=(2,1)))
init_vel = np.asmatrix(np.zeros(shape=(2,1)))
joint_pos = np.array([np.asmatrix(np.zeros(shape=(2, 1))) for i in range(n+1)])
joint_vel = np.array([np.asmatrix(np.zeros(shape=(2, 1))) for i in range(n+1)])
joint_acl = np.array([np.asmatrix(np.zeros(shape=(2, 1))) for i in range(n+1)])

# Inertial and force matrices/vectors
moi = np.zeros(n)
M = np.asmatrix(np.ones(shape=(n, n)))
M_0 = np.asmatrix(np.zeros(shape=(n, 2)))
m_0 =  mass_total * np.asmatrix(np.diag(np.ones(2)))
m = np.asmatrix(np.zeros(shape=(2, n)))
f_f = np.asmatrix(np.ones(shape=(2, 1)))
f_0 = np.asmatrix(np.ones(shape=(2, 1)))
mbar = np.zeros(n)
for i in range(n):
  # Individual link moment of inertias, modeled as solid cylinders
  moi[i] = masses[i] * (3*link_radius**2 + lengths[i]**2) / 12
  mbar[i] = masses[i] * Glengths[i] + lengths[i] * sum(masses[(i+1): n])


# Phi and Theta pos, vel, and acc vectors
Phi = np.asmatrix(np.zeros(shape=(n, 1)))
Phidot = np.asmatrix(np.zeros(shape=(n, 1)))
Phidotdot = np.asmatrix(np.zeros(shape=(n, 1)))
Theta = np.asmatrix(np.zeros(shape=(n-1, 1)))
Thetadot = np.asmatrix(np.zeros(shape=(n-1, 1)))
Thetadotdot = np.asmatrix(np.zeros(shape=(n-1, 1)))

# Friction forces
ffx = np.zeros(n)
ffy = np.zeros(n)
fft = np.zeros(n)
ffn = np.zeros(n)
sign_t = np.zeros(n)
sign_n = np.zeros(n)
u_t = np.array([np.zeros(2) for i in range(n)])
u_n = np.array([np.zeros(2) for i in range(n)])


# Plot arrays
x_t = np.array( [np.zeros(3*N_t) for i in range(n+1) ] )
y_t = np.array( [np.zeros(3*N_t) for i in range(n+1) ] )
x_i = np.zeros(n+1)
y_i = np.zeros(n+1)

# Start time cycle
for count in range(3*N_t):

  t = count * tstep
  s = vel_s * t
  Phi[0] = alpha_0 * np.cos(Bl * s) + phi_offset
  Phidot[0] = -1*alpha_0*Bl* vel_s * np.sin(Bl * s)
  Phidotdot[0] = -1*alpha_0*Bl*Bl*vel_s**2 * np.cos(Bl * s)

  for i in range(1, n):
    Theta[i-1] = -2*alpha_0*np.sin(Bn/2) * np.sin(Bl*s + Bn*i) + heading
    Thetadot[i-1] = -2*alpha_0*Bl*np.sin(Bn/2) * np.cos(Bl*s + Bn*i) * vel_s
    Thetadotdot[i-1] = 2*alpha_0*Bl*Bl*np.sin(Bn/2) * np.sin(Bl*s + Bn*i) * vel_s**2

  Phi = E * Theta + e * Phi[0]
  Phidot = E * Thetadot + e * Phidot[0]
  Phidotdot = E * Thetadotdot + e * Phidotdot[0]

  if t == 0:
    # Calculate initial joint positions
    joint_vel[0] = np.asmatrix([[np.cos(float(Phi[0]))], [np.sin(float(Phi[0]))]])

    for i in range(n+1):
      a = np.sum([lengths[j]*np.cos(Phi[j]) for j in range(i)])
      b = np.sum([lengths[j]*np.sin(Phi[j]) for j in range(i)])
      joint_pos[i] = init_pos + np.asmatrix([[a], [b]])
      x_t[i, count] = x_i[i] = float(joint_pos[i, 0])
      y_t[i, count] = y_i[i] = float(joint_pos[i, 1])

      a = np.sum([-1*lengths[j]*np.sin(Phi[j])*Phidot[j] for j in range(i)])
      b = np.sum([lengths[j]*np.cos(Phi[j])*Phidot[j] for j in range(i)])
      joint_vel[i] = joint_vel[0] + np.asmatrix([[a], [b]])
        

  else:
    
    # Update joint position
    for i in range(n+1):
      joint_pos[i] = joint_pos[i] + joint_vel[i] * tstep
      joint_vel[i] = joint_vel[i] + joint_acl[i] * tstep
      x_t[i, count] = float(joint_pos[i, 0])
      y_t[i, count] = float(joint_pos[i, 1])
  

  # Now update joint_acl

  # Calculate friction forces
  # Assumes friction point is at joint
  for i in range(n):
    u_t[i] = np.array([float(np.cos(Phi[i])), float(np.sin(Phi[i]))])
    u_n[i] = np.array([float(-1*np.sin(Phi[i])), float(np.cos(Phi[i]))])
    sign_t[i] = int(np.sign(np.inner(joint_vel[i].reshape(2), u_t[i])))
    sign_n[i] = int(np.sign(np.inner(joint_vel[i].reshape(2), u_n[i])))

    fft[i] = - mu_t * masses[i] * 9.8 * np.cos(incline) * sign_t[i]
    ffn[i] = - mu_n * masses[i] * 9.8 * np.cos(incline) * sign_n[i]

    ffx[i] = fft[i] * np.cos(Phi[i]) - ffn[i] * np.sin(Phi[i])
    ffy[i] = fft[i] * np.sin(Phi[i]) + ffn[i] * np.cos(Phi[i])


  f_f = np.asmatrix([[np.sum(ffx)], [np.sum(ffy)]])
  f_0 = np.asmatrix([[-np.sum([mbar[i]*np.cos(Phi[i])*Phidot[i]**2 for i in range(n)])],
                     [-np.sum([mbar[i]*np.sin(Phi[i])*Phidot[i]**2 for i in range(n)])]])

  # print '   t=', t, ': '
  # print 'phi: ', Phi[0]
  # print 'pos: ', joint_pos[0]
  # print 'vel: ', joint_vel[0]
  # print u_t[0], u_n[0]
  # print sign_t[0], sign_n[0]
  # print fft[0], ffn[0]
  # print ffx[0], ffy[0]
  
  for i in range(n):
    M_0[i, 0] = -1*mbar[i] * np.sin(Phi[i])
    M_0[i, 1] = mbar[i] * np.cos(Phi[i])

    m[0, i] = M_0[i, 0]
    m[1, i] = M_0[i, 1]

    for j in range(n):
      if j < i:
        M[i, j] = lengths[i] * mbar[j] * np.cos(Phi[j] - Phi[i])
      elif i == j:
        M[i, j] = masses[i] * Glengths[i]**2 + moi[i] + lengths[i]**2 * np.sum(masses[(i+1):])
      elif j > i:
        M[i, j] = mbar[i] * lengths[j] * np.cos(Phi[j] - Phi[i])


  # print "m: "
  # print m, '\n'
  # print "m_0: "
  # print m_0, '\n'
  # print "f0"
  # print f_0
  # print
  # print "ff"
  # print f_f
  # print
  # print "phidotdot"
  # print Phidotdot
  # print

#  print t, joint_pos[0].reshape(2)

  # Calculate joint acceleration from forces.
#  print '  t=', t
  # print m
  # print Phidotdot
  # print np.linalg.inv(m_0)
  # joint_acl[0] = -np.linalg.inv(m_0) * m * Phidotdot - np.linalg.inv(m_0) * (f_0 + f_f) - g
  # print joint_acl[0]
  # print np.linalg.inv(m_0) * m * Phidotdot
  # print np.linalg.inv(m_0) * (f_0 + f_f)m * Phidotdot - f_f - f_0
  a = - np.sum([mbar[i]*np.sin(Phi[i])*Phidotdot[i] for i in range(n)])
  b = np.sum([mbar[i]*np.cos(Phi[i])*Phidotdot[i] for i in range(n)])
#  print np.asmatrix([[a], [b]]) * -(1/mass_total)
  joint_acl[0] = -(1/mass_total) * (np.asmatrix([[a], [b]]) - f_f - f_0) - g
  # print (np.asmatrix([[a], [b]]) - f_f + f_0)
  # print joint_acl[0], '\n'
  # calc other joint acl 
  for i in range(n+1):
    a = - np.sum([lengths[j]*(np.sin(Phi[j])*Phidotdot[j] +
                              np.cos(Phi[j]) * Phidot[j]**2) for j in range(i)])
    b = np.sum([lengths[j]*(np.cos(Phi[j])*Phidotdot[j] -
                              np.sin(Phi[j]) * Phidot[j]**2) for j in range(i)])
    joint_acl[i] = joint_acl[0] + np.asmatrix([[a], [b]])

  # print 'p0_acl: ', joint_acl[0]
  # print

  # End of time cycle

plt.plot(x_i, y_i)
plt.plot(x_t[0], y_t[0])
plt.plot(x_t[n], y_t[n])
plt.ylabel('y, meters')
plt.xlabel('x, mmeters')
plt.axis('equal')
plt.grid()
plt.show()
  

























# for t in range(0, period * 2, tstep):

  # s(t)

'''
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


'''


# How to plot
# plt.plot(x_i, y_i)
# plt.ylabel('y, meters')
# plt.xlabel('x, mmeters')
# plt.axis('equal')
# plt.grid()
# plt.show()
