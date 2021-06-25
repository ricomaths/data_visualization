from aux_functions import *

# # input data
Pe=0.01
n_r=1000
n_theta=750

# initializing
r = np.zeros([n_r*n_theta,1],dtype=np.float64)
theta_angle=np.zeros([n_r*n_theta,1],dtype=np.float64)
theta_in = np.zeros([n_r*n_theta,1],dtype=complex)
theta_out = np.zeros([n_r*n_theta,1],dtype=complex)

# t=100
for t in np.arange(50,1050,50):
    # # inner region
    r_max_inner=10
    [r[:,0],theta_angle[:,0],theta_in[:,0]] = inner_region(Pe,t,r_max_inner,n_r,n_theta)
    temperature_2D_polar_2vtk(r,theta_angle,n_r,n_theta,theta_in,"inner_region_t"+str(t))

    # # outer region
    r_max_outer=300
    [r[:,0],theta_angle[:,0],theta_out[:,0]] = outer_region(Pe,t,r_max_outer,n_r,n_theta)
    temperature_2D_polar_2vtk(r,theta_angle,n_r,n_theta,theta_out,"outer_region_t"+str(t))