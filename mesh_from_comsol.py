import numpy as np

def reading_mesh_comsol(path):
    with open(path,'r') as f:
        lines_list = f.readlines()
        n_points = int(lines_list[17].strip().split()[0])
        points_cloud = np.zeros([n_points,2],dtype=float)
        for i in np.arange(0,n_points,1):
            points_cloud[i,0] = lines_list[21+i].strip().split()[0]
            points_cloud[i,1] = lines_list[21+i].strip().split()[1]
        n_elements = int(lines_list[21+n_points+9].strip().split()[0])
        connectivity_matrix = np.zeros([n_elements,3],dtype=int)
        for i in np.arange(0,n_elements,1):
            for j in np.arange(0,3,1):
                connectivity_matrix[i,j] = lines_list[21+n_points+11+i].strip().split()[j]
    return n_points,points_cloud,n_elements,connectivity_matrix
