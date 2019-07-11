import numpy as np
import math
import matplotlib.pyplot as plt

# robot parameters
mu_t = 0.01
mu_n = 0.9
n = 8
link_mass = 1.
link_length = 1.
link_radius = 0.1
lG = 0.5 # distance to center of gravity of link
lf = 0.1 # distance to friction point (wheels)
masses = np.repeat(link_mass, n)
lengths = np.repeat(link_length, n)
Glengths = np.repeat(lG, n)
flengths = np.repeat(lf, n)
L = np.sum(lengths)
mass_total = np.sum(masses)

# motion parameters
incline = 0 #(np.pi/180) * 10
alpha_0 = np.pi / 6.
Kn = 1
vel_s = 1.2
phi_offset = 0
heading = 0
Bl = 2*Kn*np.pi/L
Bn = 2*Kn*np.pi/n
period = (L / vel_s)
N_t = 100
tstep = period / N_t
g = np.asmatrix([[9.8*np.sin(incline)], [0]])
sim_time = 3*N_t

  
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
center_pos = np.asmatrix(np.zeros(shape=(2, 1)))
center_vel = np.asmatrix(np.zeros(shape=(2, 1)))


# Inertial and force, torque matrices/vectors
moi = np.zeros(n)
M = np.asmatrix(np.ones(shape=(n, n)))
M_0 = np.asmatrix(np.zeros(shape=(n, 2)))
m_0 =  mass_total * np.asmatrix(np.diag(np.ones(2)))
m = np.asmatrix(np.zeros(shape=(2, n)))
f_f = np.asmatrix(np.ones(shape=(2, 1)))
f_0 = np.asmatrix(np.ones(shape=(2, 1)))
T = np.asmatrix(np.zeros(shape=(n-1, 1)))
T_f = np.asmatrix(np.zeros(shape=(n, 1)))
T_0 = np.asmatrix(np.zeros(shape=(n, 1)))
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
x_t = np.array( [np.zeros(sim_time) for i in range(n+1) ] )
y_t = np.array( [np.zeros(sim_time) for i in range(n+1) ] )
x_i = np.zeros(n+1)
y_i = np.zeros(n+1)
v_t = np.zeros(sim_time)
c_t = np.array( [np.zeros(sim_time) for i in range(2)] )
T_t = np.array( [np.zeros(sim_time) for i in range(n-1)] )
time = np.zeros(3*N_t)

# Start time cycle
for count in range(sim_time):

  t = count * tstep
  time[count] = t
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

  rotvec = np.array( [np.asmatrix([[np.cos(float(Phi[i]))], [np.sin(float(Phi[i]))]])
                      for i in range(n)] )

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


  # average link angle
  g_pos = np.array( [ (joint_pos[i] + Glengths[i] * rotvec[i]) for i in range(n) ] )
  phibar = float(sum(Phi)) / n
  c_t[0, count] = sum([masses[i]*g_pos[i, 0] for i in range(n)]) / mass_total
  c_t[1, count] = sum([masses[i]*g_pos[i, 1] for i in range(n)]) / mass_total
  c_vel = sum( [masses[i]*joint_vel[i] for i in range(n)] ) / mass_total
  v_t[count] = np.inner(c_vel.reshape(2), np.array([np.cos(phibar), np.sin(phibar)]))
  
  
  # Now update joint_acl

  # Calculate friction and torque forces
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

    T_f[i] = (np.sin(Phi[i])*(lengths[i]*np.sum(ffx[i+1:])+flengths[i]*ffx[i]) -
              np.cos(Phi[i])*(lengths[i]*np.sum(ffy[i+1:])+flengths[i]*ffy[i]))
    T_0[i] = (lengths[i]*sum([mbar[k]*np.sin(Phi[k]-Phi[i])*Phidot[k]**2 for k in range(i+1, n)]) -
              mbar[i]*sum([lengths[k]*np.sin(Phi[k]-Phi[i])*Phidot[k]**2 for k in range(i)]))


  f_f = np.asmatrix([[np.sum(ffx)], [np.sum(ffy)]])
  f_0 = np.asmatrix([[-np.sum([mbar[i]*np.cos(Phi[i])*Phidot[i]**2 for i in range(n)])],
                     [-np.sum([mbar[i]*np.sin(Phi[i])*Phidot[i]**2 for i in range(n)])]])
  
  
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


    # Calculate joint acceleration from forces.
  joint_acl[0] = np.linalg.inv(m_0) * m * Phidotdot + np.linalg.inv(m_0) * (f_0 + f_f) - g

  DT = T_0 + T_f + M_0*(joint_acl[0] + g) + M*Phidotdot

  # Calculate joint torques
  for i in range(n-1):
    if(i == 0):
      T[i] = - DT[i]
    else:
      T[i] = T[i-1] - DT[i]
    T_t[i, count] = T[i]      


  # calc other joint acl 
  for i in range(n+1):
    a = - np.sum([lengths[j]*(np.sin(Phi[j])*Phidotdot[j] +
                              np.cos(Phi[j]) * Phidot[j]**2) for j in range(i)])
    b = np.sum([lengths[j]*(np.cos(Phi[j])*Phidotdot[j] -
                              np.sin(Phi[j]) * Phidot[j]**2) for j in range(i)])
    joint_acl[i] = joint_acl[0] + np.asmatrix([[a], [b]])


    
    
  # print '   t=', t, ': '
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
  # print 'phi: '
  # print Phi[0]
  # print 'pos: '
  # print joint_pos[0]
  # print 'vel: '
  # print joint_vel[0]
  # print 'unit vecs: '
  # print u_t[0], u_n[0]
  # print 'signs: '
  # print sign_t[0], sign_n[0]
  # print 'fft, ffn: '
  # print fft[0], ffn[0]
  # print 'ffx, ffy'
  # print ffx[0], ffy[0]
  # print 'f_f'
  # print f_f
  # print 'f_0'
  # print f_0
  # print 'p0_acl: ', joint_acl[0]
  # print
  # print


  # End of time cycle

v_avg = sum(v_t) / sim_time
print 'avg vel_t: ', v_avg
  
plt.plot(x_i, y_i)
plt.plot(x_t[0], y_t[0])
plt.plot(x_t[n], y_t[n])
plt.plot(c_t[0], c_t[1])
plt.ylabel('y, meters')
plt.xlabel('x, mmeters')
plt.axis('equal')
plt.grid()
plt.show()

for i in range(n-1):
  plt.plot(time, T_t[i])
plt.show()

plt.plot(time, v_t)
plt.show()
