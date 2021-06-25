import numpy as np
import scipy.special as sp

def inner_region(Pe,t,r_max,n_r,n_theta):
    # # initializing
    r_inner = np.zeros(n_r*n_theta)
    theta_angle = np.zeros(n_r*n_theta)
    # # meshing
    for i in np.arange(0,n_r,1):
        for j in np.arange(0,n_theta,1):
            theta_angle[i*n_theta+j] = 2*np.pi*j/n_theta
            r_inner[i*n_theta+j] = 1+(r_max-1)*i/n_r

    # # boundary condition
    q= 2*np.pi
    q_m1 = complex(1,1)/2.0
    q_1 = complex(1,-1)/2.0

    # # inner solution
    theta_inn = -q/(2*np.pi) *np.log(r_inner) -q/(4*np.pi) *(sp.exp1(Pe**2/4*t)+np.log(Pe**2/16)+2*np.euler_gamma)\
        +q_m1*np.exp(-theta_angle*complex(0,1))/r_inner + q_1*np.exp(theta_angle*complex(0,1))/r_inner
    return r_inner,theta_angle,theta_inn

def temperature_2D_polar_2vtk(r,theta_angle,n_r,n_theta,image,name_output):
    # # creating the vtk file
    vtk_file=open(name_output + '.vtk','w')
    # # Identifier
    vtk_file.write("# vtk DataFile Version 3.0\n")
    # # Header
    vtk_file.write("Temperature field\n")
    # # Format
    vtk_file.write("ASCII\n")
    # # Data set structure
    vtk_file.write("DATASET STRUCTURED_GRID\n")
    # # Mesh data
    vtk_file.write("DIMENSIONS "+str(n_r)+" "+str(n_theta)+" 1\n")
    vtk_file.write("POINTS " + str(n_r*n_theta)+" double\n")
    for i in np.arange(0,n_r*n_theta,1):
        vtk_file.write(str(r[i,0]*np.cos(theta_angle[i,0]))+" "+str(r[i,0]*np.sin(theta_angle[i,0]))+" 0\n")
    # # Real Temperature data
    vtk_file.write("POINT_DATA "+str(n_r*n_theta)+"\n")
    vtk_file.write("SCALARS "+ "Modulus_of_the_Temperature "+"DOUBLE\n")
    vtk_file.write("LOOKUP_TABLE default\n")
    for i in np.arange(0,n_r*n_theta,1):
        vtk_file.write(str(abs(image[i,0]))+"\n")
    vtk_file.write("SCALARS "+ "Argument_of_the_Temperature "+"DOUBLE\n")
    vtk_file.write("LOOKUP_TABLE default\n")
    for i in np.arange(0,n_r*n_theta,1):
        vtk_file.write(str(np.angle(image[i,0]))+"\n")
    # # closing the file
    vtk_file.close()

# # auxiliary integrals
def integral_exp_sqrt_ab(t,Pe,r):
    a = Pe**2/4*t+r**2/(4*t)
    b = r*Pe/2
    integral = 1/(3951360*a**4)*(-20152*b**2*np.exp(-np.sqrt(126 - 7*np.sqrt(30))*b/14)*a**4*(np.sqrt(30) + 24370/2519)*sp.exp1( a - np.sqrt(126 - 7*np.sqrt(30))*b/14) - 20152*b**2*np.exp(np.sqrt(126 - 7*np.sqrt(30))*b/14)*a**4*(np.sqrt(30) + 24370/2519)*sp.exp1( a + np.sqrt(126 - 7*np.sqrt(30))*b/14) + 20152*b**2*(np.sqrt(30) - 24370/2519)*a**4*np.exp(-np.sqrt(126 + 7*np.sqrt(30))*b/14)*sp.exp1( a - np.sqrt(126 + 7*np.sqrt(30))*b/14) + 20152*np.exp(np.sqrt(126 + 7*np.sqrt(30))*b/14)*b**2*(np.sqrt(30) - 24370/2519)*a**4*sp.exp1( a + np.sqrt(126 + 7*np.sqrt(30))*b/14) + 20152*a**4*b**2*(np.sqrt(30) + 24370/2519)*np.exp(-np.sqrt(126 - 7*np.sqrt(30))*b/14)*sp.exp1( b - np.sqrt(126 - 7*np.sqrt(30))*b/14) + 20152*np.exp(np.sqrt(126 - 7*np.sqrt(30))*b/14)*a**4*b**2*(np.sqrt(30) + 24370/2519)*sp.exp1( b + np.sqrt(126 - 7*np.sqrt(30))*b/14) - 20152*a**4*np.exp(-np.sqrt(126 + 7*np.sqrt(30))*b/14)*(np.sqrt(30) - 24370/2519)*b**2*sp.exp1( b - np.sqrt(126 + 7*np.sqrt(30))*b/14) - 20152*np.exp(np.sqrt(126 + 7*np.sqrt(30))*b/14)*a**4*(np.sqrt(30) - 24370/2519)*b**2*sp.exp1( b + np.sqrt(126 + 7*np.sqrt(30))*b/14) + (3951360*a**5 + 3951360*a**4 + (245*b**6 + 39480*b**4)*a**3 + (-245*b**6 - 39480*b**4)*a**2 + 490*b**6*a - 1470*b**6)*np.exp(-a) - 245*((b**5 - b**4 + 1142/7*b**3 - 1170/7*b**2 + 16128*b + 16128)*np.exp(-b) + (b**4 + 1128/7*b**2 + 239168/49)*(sp.exp1( a) - sp.exp1( b))*b**2)*a**4)
    return integral

def integral_exp_z_sqrt_ab(t,Pe,r):
    a = Pe**2/4*t+r**2/(4*t)
    b = r*Pe/2
    integral = 1/(np.sqrt(126+7*np.sqrt(30))*np.sqrt(126-7*np.sqrt(30))*a**3)*(-(107*a**3*np.sqrt(126 + 7*np.sqrt(30))*(np.sqrt(30) + 3705/214)*np.exp(-np.sqrt(126 - 7*np.sqrt(30))*b/14)*b**3*sp.exp1( a - np.sqrt(126 - 7*np.sqrt(30))*b/14))/5040 + (107*a**3*np.exp(np.sqrt(126 - 7*np.sqrt(30))*b/14)*np.sqrt(126 + 7*np.sqrt(30))*(np.sqrt(30) + 3705/214)*b**3*sp.exp1( a + np.sqrt(126 - 7*np.sqrt(30))*b/14))/5040 + (107*b**3*np.sqrt(126 - 7*np.sqrt(30))*np.exp(-np.sqrt(126 + 7*np.sqrt(30))*b/14)*(np.sqrt(30) - 3705/214)*a**3*sp.exp1( a - np.sqrt(126 + 7*np.sqrt(30))*b/14))/5040 - (107*np.exp(np.sqrt(126 + 7*np.sqrt(30))*b/14)*b**3*np.sqrt(126 - 7*np.sqrt(30))*(np.sqrt(30) - 3705/214)*a**3*sp.exp1( a + np.sqrt(126 + 7*np.sqrt(30))*b/14))/5040 + (107*b**3*np.sqrt(126 + 7*np.sqrt(30))*np.exp(-np.sqrt(126 - 7*np.sqrt(30))*b/14)*(np.sqrt(30) + 3705/214)*a**3*sp.exp1( b - np.sqrt(126 - 7*np.sqrt(30))*b/14))/5040 - (107*np.exp(np.sqrt(126 - 7*np.sqrt(30))*b/14)*b**3*np.sqrt(126 + 7*np.sqrt(30))*(np.sqrt(30) + 3705/214)*a**3*sp.exp1( b + np.sqrt(126 - 7*np.sqrt(30))*b/14))/5040 + np.sqrt(126 - 7*np.sqrt(30))*(-(107*b**3*np.exp(-np.sqrt(126 + 7*np.sqrt(30))*b/14)*(np.sqrt(30) - 3705/214)*a**3*sp.exp1( b - np.sqrt(126 + 7*np.sqrt(30))*b/14))/5040 + (107*np.exp(np.sqrt(126 + 7*np.sqrt(30))*b/14)*b**3*(np.sqrt(30) - 3705/214)*a**3*sp.exp1( b + np.sqrt(126 + 7*np.sqrt(30))*b/14))/5040 + np.sqrt(126 + 7*np.sqrt(30))*((a**5 + 2*a**4 + (-b**2/2 + 2)*a**3 + (-1/4032*b**6 - 47/2352*b**4)*a**2 + b**6*a/4032 - b**6/2016)*np.exp(-a) + (((b**5 - b**4 + 578/7*b**3 - 2016*b**2 - 8064*b - 8064)*np.exp(-b) + (sp.exp1( a) - sp.exp1( b))*b**4*(b**2 + 564/7))*a**3)/4032)))
    return integral

def integral_1(t,Pe,r):
    a = Pe**2/4*t+r**2/(4*t)
    b = r*Pe/2
    integral = 1/b**2*(2*integral_exp_sqrt_ab(t,Pe,r)+np.exp(-a)*a*np.sqrt(a**2-b**2)-integral_exp_z_sqrt_ab(t,Pe,r))
    return integral

def integral_exp_sqrt_binfinity(Pe,r):
    b = r*Pe/2
    integral=-(2519*b**2*np.exp(-np.sqrt(126 - 7*np.sqrt(30))*b/14)*(np.sqrt(30) + 24370/2519)*sp.exp1( b - np.sqrt(126 - 7*np.sqrt(30))*b/14))/493920 - (2519*np.exp(np.sqrt(126 - 7*np.sqrt(30))*b/14)*b**2*(np.sqrt(30) + 24370/2519)*sp.exp1( b + np.sqrt(126 - 7*np.sqrt(30))*b/14))/493920 + (2519*b**2*(np.sqrt(30) - 24370/2519)*np.exp(-np.sqrt(126 + 7*np.sqrt(30))*b/14)*sp.exp1( b - np.sqrt(126 + 7*np.sqrt(30))*b/14))/493920 + (2519*b**2*(np.sqrt(30) - 24370/2519)*np.exp(np.sqrt(126 + 7*np.sqrt(30))*b/14)*sp.exp1( b + np.sqrt(126 + 7*np.sqrt(30))*b/14))/493920 + ((245*b**5 - 245*b**4 + 39970*b**3 - 40950*b**2 + 3951360*b + 3951360)*np.exp(-b))/3951360 - b**2*(b**4 + 1128/7*b**2 + 239168/49)*sp.exp1( b)/16128
    return integral

def integral_exp_z_sqrt_binfinity(Pe,r):
    b = r*Pe/2
    integral=1/(4032*np.sqrt(126-7*np.sqrt(30))*np.sqrt(126+7*np.sqrt(30)))*(-(428*np.sqrt(126 + 7*np.sqrt(30))*b**3*np.exp(-np.sqrt(126 - 7*np.sqrt(30))*b/14)*(np.sqrt(30) + 3705/214)*sp.exp1( b - np.sqrt(126 - 7*np.sqrt(30))*b/14))/5 + (428*np.exp(np.sqrt(126 - 7*np.sqrt(30))*b/14)*np.sqrt(126 + 7*np.sqrt(30))*b**3*(np.sqrt(30) + 3705/214)*sp.exp1( b + np.sqrt(126 - 7*np.sqrt(30))*b/14))/5 + np.sqrt(126 - 7*np.sqrt(30))*(428*(np.sqrt(30) - 3705/214)*b**3*np.exp(-np.sqrt(126 + 7*np.sqrt(30))*b/14)*sp.exp1( b - np.sqrt(126 + 7*np.sqrt(30))*b/14)/5 - (428*np.exp(np.sqrt(126 + 7*np.sqrt(30))*b/14)*(np.sqrt(30) - 3705/214)*b**3*sp.exp1( b + np.sqrt(126 + 7*np.sqrt(30))*b/14))/5 + np.sqrt(126 + 7*np.sqrt(30))*((-b**5 + b**4 - 578/7*b**3 + 2016*b**2 + 8064*b + 8064)*np.exp(-b) + sp.exp1( b)*b**4*(b**2 + 564/7))))
    return integral

def integral_2(Pe,r):
    b=r*Pe/2
    integral = 1/b**2*(integral_exp_z_sqrt_binfinity(Pe,r)-2*integral_exp_sqrt_binfinity(Pe,r))
    return integral

# base integral (case n=0)
def integral_base(t,Pe,r):
    integral = integral_1(t,Pe,r) + integral_2(Pe,r)
    return integral

# integral without numerator
def integral_no_numerator(t,Pe,r):
    a = Pe**2/4*t+r**2/(4*t)
    b = r*Pe/2
    integral = 2/Pe**2*np.exp(-a) - 2/Pe**2*(-np.exp(-a)*np.sqrt(a**2-b**2)+integral_exp_sqrt_ab(t,Pe,r)) + 2/Pe**2*integral_exp_sqrt_binfinity(Pe,r)
    return integral

# integral with numerator tau**2
def integral_numerator_tau2(t,Pe,r):
    integral = 4/r**2*(Pe**2/4*integral_no_numerator(t,Pe,r)-np.exp(-Pe**2/4*t-r**2/(4*t)))
    return integral

def integral_outer_tinfinity(t,Pe,r,n_max):
    integral = np.zeros([len(r),n_max])
    integral[:,0] = integral_base(t,Pe,r)
    integral[:,1] = integral_numerator_tau2(t,Pe,r)
    for i in np.arange(2,n_max,1):
        integral[:,i] = Pe**2/r**2*(integral[:,i-2] + 4*(i-1)/Pe**2*integral[:,i-1]-4/Pe**2*np.exp(-Pe**2*t/4-r**2/(4*t))/t**(i-1))
    return integral

def integral_outer_0t(t,Pe,r,n_max):
    z = r*Pe/2
    integral = np.zeros([len(r),n_max])
    for i in np.arange(0,n_max,1):
        integral[:,i] =  2*(2/z)**i*(Pe**2/4)**i*sp.kn(i,z)-integral_outer_tinfinity(t,Pe,r,n_max)[:,i]
    return integral 

def outer_region(Pe,t,r_max_outer,n_r,n_theta):
    r_min=1/np.sqrt(Pe**2+1/t)

    # # # # # mesh
    r_outer = np.zeros(n_r*n_theta)
    theta_angle = np.zeros(n_r*n_theta)

    for i in np.arange(0,n_r,1):
        for j in np.arange(0,n_theta,1):
            theta_angle[i*n_theta+j] = 2*np.pi*j/n_theta
            r_outer[i*n_theta+j] = r_min+(r_max_outer-r_min)*i/n_r


    # # # boundary condition 
    q0=2*np.pi
    q_m1 = complex(1,1)/2.0
    q_1 = complex(1,-1)/2.0
    n_max=2

    integral_convolution = integral_outer_0t(t,Pe,r_outer,n_max)

    # # # outer solution

    theta_outer = np.exp(r_outer*Pe*np.cos(theta_angle)/2) *(q0/(4*np.pi)*integral_convolution[:,0] +r_outer/2*1/2*integral_convolution[:,1]\
        *(q_m1*np.exp(-theta_angle*complex(0,1)) + q_1*np.exp(theta_angle*complex(0,1)) )  )

    return r_outer,theta_angle,theta_outer