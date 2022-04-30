# name: multipole_bc_solver
# meaning: multipole, boundary_condition solver
# usage:    using multipole method,
#           knowing the boundary condition of the vibrating body,
#           solving the coefficences of the multipoles
# inputs:   a specified frequency OMEGA
#           a sery of boundary condition: 
#               i.e.    the modal displacements in M different sample point
#                       the normal vector of the boundary of this vibratin body in this M different point 
#           a sery of locations of multipole
#               totally N different multipole (in N different location)
#               each location in Cartesean coordination (i.e. [x, y, z])
# outputs:  a sery of coefficience of this multipole:
#               totally 4*N (since every multipole is expanded with 2-degree, 4 coeffiences)

import sys
from typing import Any
import numpy as np
from pandas import array
import scipy.linag

## macro:
RHO = 1.2041    # density of air in kg/m^3, for compute the boudary condition
OMEGA = 2763    # a specified frequence, in Hz, will be a constant if the modal is specified
c = 340         # in m/s velocity of sound in the air.

class Multipole:
    def __init__(self, location, coeffiencience) -> None:
        self.location=location                              ## Cartesian coordinates [x, y, z] 
        self.coeffience=coeffiencience                      ## 2-degree multipole expansion, so we need 4 coefficience
        pass

k = OMEGA / c       # wave number, a constant for specified modal
def get_sampling_location_series(file)->Any:
    'This function should read in the modal displacement'
    'of M Boundary Conditions[BCs]'
    "这是一个用于读入“模态采样点位置”的函数"
    ... 
    filein
    ...
    sampling_location_series = ...
    return sampling_location_series
    ## this sampling_location_series should be a N * 1*3 vector?
    ## 返回的应该是有N个1*3向量的数组
    
def get_modal_displacement_series(file)->Any:
    'This function should read in the modal displacement'
    'of M Boundary Conditions[BCs]'
    "这是一个用于读入“模态边界振动情况”的函数"
    ... ## file_in operation
    filein
    ... ## read operation
    modal_displacement_series = ... 
    return modal_displacement_series 
    ## this modal_displacement_series should be a M * 1*3 vector?
    ## 返回的应该是有M个1*3向量的数组
    
def get_boundary_normal_vector_series(file)->Any:
    'This function should read in the normal vector of the boundary of the body'
    "这是一个读入物体轮廓法向的函数"
    ...
    filein
    ...
    normal_vector_series = ... 
    return normal_vector_series
    ## the normal_vecotr_series should be a M * 1 * 3 vector
    ## 返回的应该是有M个1*3向量的数组

def  get_multipole_location_series(multipole_location_file)->Any:
    'This function should read in the location of N multipole'
    "这是一个读入多级源位置的函数"
    ...
    multipole_location_file
    ...
    multipole_location_series = ...
    return multipole_location_series
    ## the multipole_location_series should be a N * 1 * 3 vector
    ## 返回的应该是有N个1*3向量的数组
    
def construct_rhs(normal_vector_series, modal_displacement_series)-> Any:
    "construct right hand side"
    "This function takes in normal_vector and modal_displacement"
    "gives out the boundary condition"
    '构造方程右边，读入M * 1 * 3的“法向向量”和M * 1 *3的“模态振动？？”，输出M*1的边界条件bc'
    if len(normal_vector_series) == len(modal_displacement_series):
        M = len(normal_vector_series)
    else:
        print("the size of normal_vecor and the modal_vector do not match, construct_rhs error")
        return
    a = np.array(normal_vector_series)
    b = np.array(modal_displacement_series)
    c = (a * b).sum(axis = 1)
    bc = OMEGA ** 2 * RHO * c           ## bc means boundary_condition
    if np.shape(bc)==(M,):
        return bc
    else:
        print("the return size of bc is not M*1, construct_rhs error")
    
def normal_derivative_S_0_0(k, sampling_loc, multipole_loc, normal_vector)->Any:
    "takes in one sampling location[1*3],one multipole location[1*3], one normal vector[1*3] and k"
    "return the normal derivative of S_0_0 in ONLY ONE POINT"
    "MUST note that this normal_vector is specified to only one point, not M point as in the funct construct_rhs"
    '读入一个采样点位置,一个多级源位置，一个采样点位置的法向，k，计算这个S_0^0的法向导数【梯度和法向的内积】'
    r = np.sqrt(np.sum((np.array(sampling_loc) - np.array(multipole_loc))**2))          # r is the dist from sample location to multipole location
    delta = np.dot(np.array(sampling_loc) - np.array(multipole_loc), normal_vector)     # delta0 is 
    res = -1/2 * np.sqrt(1/np.pi) * (1/(k * r**3)) * (0.0 + 1.0j) * np.exp(0.0 - 1.0j * k * r) * (1.0 + 1.0j * k * r) * delta
    return res

def  normal_derivative_S_1_minus1(k, sampling_loc, multipole_loc, normal_vector)->Any:
    "takes in one sampling location[1*3],one multipole location[1*3], one normal vector[1*3] and k"
    "return the normal derivative of S_1_-1 in ONLY ONE POINT"
    "MUST note that this normal_vector is specified to only one point, not M point as in the funct construct_rhs"
    '读入一个采样点位置,一个多级源位置，一个采样点位置的法向，k，计算这个S_1^-1的法向导数【梯度和法向的内积】'
    r = np.sqrt(np.sum((np.array(sampling_loc) - np.array(multipole_loc))**2))                          # r is the dist from sample location to multipole location
    delta0 = np.dot(np.array(sampling_loc) - np.array(multipole_loc), normal_vector)                    # delta0 
    diff = np.array(sampling_loc) - np.array(multipole_loc)                                             # 'diff' means difference vector
    diff_x = diff[0]                                                                                    # x_i - x_{0j}
    diff_y = diff[1]                                                                                    # y_i - y_{0j}
    diff_z = diff[2]                                                                                    # z_i - z_{0j}
    item = [r ** 2 - diff_x ** 2 + 1.0j * diff_x * diff_y, -(diff_y * diff_x + 1.0j * (r**2 - diff_y**2)), diff_z * (diff_x - 1.0j * diff_y)]
    delta1 = np.dot(np.array(normal_vector, item))
    #the var 'res' below can only looking for the formula in markdown to check and debug...
    res = -1/2 * np.sqrt(3 / (2 * np.pi)) * np.exp(0.0 - 1.0j*k*r) * (1/(k**2*r**5)) * ((k * r - 1.0j)*delta1 + (-2*k*r + 1.0j * (2 - k**2*r**2))*(diff_x-1.0j*diff_y)*delta0)   
    return res

def  normal_derivative_S_1_0(k, sampling_loc, multipole_loc, normal_vector)->Any:
    "takes in one sampling location[1*3],one multipole location[1*3], one normal vector[1*3] and k"
    "return the normal derivative of S_1_0 in ONLY ONE POINT"
    "MUST note that this normal_vector is specified to only one point, not M point as in the funct construct_rhs"
    '读入一个采样点位置,一个多级源位置，一个采样点位置的法向，k，计算这个S_1^0的法向导数【梯度和法向的内积】'
    r = np.sqrt(np.sum((np.array(sampling_loc) - np.array(multipole_loc))**2))                          # r is the dist from sample location to multipole location
    delta0 = np.dot(np.array(sampling_loc) - np.array(multipole_loc), normal_vector)                    # delta0 
    diff = np.array(sampling_loc) - np.array(multipole_loc)                                             # 'diff' means difference vector
    diff_x = diff[0]                                                                                    # x_i - x_{0j}
    diff_y = diff[1]                                                                                    # y_i - y_{0j}
    diff_z = diff[2]                                                                                    # z_i - z_{0j}
    item = [diff_z * diff_x, diff_z * diff_y, r**2 - diff_z**2]
    delta1 = np.dot(np.array(normal_vector, item))                                                      # delta1
    #the var 'res' below can only looking for the formula in markdown to check and debug...
    res = -1/2 * np.sqrt(3 / (2 * np.pi)) * np.exp(0.0 - 1.0j*k*r) * (1/(k**2*r**5)) * ((k * r - 1.0j)*delta1 + (-2*k*r + 1.0j * (2 - k**2*r**2))*(diff_z)*delta0)   
    return res

def  normal_derivative_S_1_1(k, sampling_loc, multipole_loc, normal_vector)->Any:
    "takes in one sampling location[1*3],one multipole location[1*3], one normal vector[1*3] and k"
    "return the normal derivative of S_1_1 in ONLY ONE POINT"
    "MUST note that this normal_vector is specified to only one point, not M point as in the funct construct_rhs"
    '读入一个采样点位置,一个多级源位置，一个采样点位置的法向，k，计算这个S_1^1的法向导数【梯度和法向的内积】'
    r = np.sqrt(np.sum((np.array(sampling_loc) - np.array(multipole_loc))**2))                          # r is the dist from sample location to multipole location
    delta0 = np.dot(np.array(sampling_loc) - np.array(multipole_loc), normal_vector)                    # delta0 
    diff = np.array(sampling_loc) - np.array(multipole_loc)                                             # 'diff' means difference vector
    diff_x = diff[0]                                                                                    # x_i - x_{0j}
    diff_y = diff[1]                                                                                    # y_i - y_{0j}
    diff_z = diff[2]                                                                                    # z_i - z_{0j}
    item = [r ** 2 - diff_x ** 2 - 1.0j * diff_x * diff_y, -diff_y * diff_x + 1.0j * (r**2 - diff_y**2), diff_z * (diff_x + 1.0j * diff_y)]
    delta1 = np.dot(np.array(normal_vector, item))
    #the var 'res' below can only looking for the formula in markdown to check and debug...
    res = -1/2 * np.sqrt(3 / (2 * np.pi)) * np.exp(0.0 - 1.0j*k*r) * (1/(k**2*r**5)) * ((k * r - 1.0j)*delta1 + (-2*k*r + 1.0j * (2 - k**2*r**2))*(diff_x+1.0j*diff_y)*delta0)   
    return res

def construct_coefficiency_matrix(k, sampling_location_series, multipole_location_series, normal_vector_series)->Any:
    if len(sampling_location_series) == len(normal_vector_series):
        M = len(sampling_location_series)
    else:
        print("The size of sampling_location_series and normal_vector_series does not match, consturct_coefficiency_matrix error")
        return
    N = len(multipole_location_series)
    A = np.zeros((M, 4 * N), dtype=complex)
    for i in range(0, M):                                                   # 不需要M+1, i in range(0, M)是从0遍历到M-1,正好M行 
        for j in range(0, 4 * N):
            sampling_loc_i = sampling_location_series[i]
            normal_vec_i = normal_vector_series[i]
            multipole_loc_j = multipole_location_series[j]
            if j % 4 == 0:
                A[i][j] = normal_derivative_S_0_0(k=k,sampling_loc=sampling_loc_i, multipole_loc=multipole_loc_j, normal_vector=normal_vec_i)
            if j % 4 == 1:
                A[i][j] = normal_derivative_S_1_minus1(k=k, sampling_loc=sampling_loc_i, multipole_loc=multipole_loc_j, normal_vector=normal_vec_i)
            if j % 4 == 2:
                A[i][j] = normal_derivative_S_1_0(k=k, sampling_loc=sampling_loc_i, multipole_loc=multipole_loc_j, normal_vector=normal_vec_i)
            if j % 4 == 3:
                A[i][j] = normal_derivative_S_1_1(k=k, sampling_loc=sampling_loc_i, multipole_loc=multipole_loc_j, normal_vector=normal_vec_i)

    return A

def solver(file, file2, opt):                                            # if multipole file, try like   solver(Afile, Bfile, ..., file2)
    if opt >= 3 or opt <= 0:
        print("unsupported solution option, solver error")
    sampling_loc_sr = get_sampling_location_series(file=file)       # if multipole file, try like   get_sampling_location_series(file=Afile), 
    modal_disp_sr = get_modal_displacement_series(file=file)        #                               get_modal_displacment_series(file=Bfile),
    boundary_n_sr = get_boundary_normal_vector_series(file=file)    #                               e.t.c
    
    multipole_loc_sr = get_multipole_location_series(file2=file2)
    
    M = len(sampling_loc_sr)
    N = len(multipole_loc_sr)
    rhs  = construct_rhs(normal_vector_series=boundary_n_sr, modal_displacement_series=modal_disp_sr)
    coefficiency_matrix = construct_coefficiency_matrix(k, sampling_location_series=sampling_loc_sr,multipole_location_series=multipole_loc_sr, normal_vector_series=boundary_n_sr)
    if opt == 1:
    #solving plan A
        c = np.linalg.solve(coefficiency_matrix, rhs)
    elif opt == 2:
    # solving plan B
        c, res, rnk, s = scipy.linalg.lstsq(coefficiency_matrix, rhs)
    # check the length of c:
    if np.shape(c) == (4 * N, ):
        ... 
        sys.fileout 
        'coefficiency_file'
        c 
        fileout 
        ... # write out the c file
        return c
    else:
        print("solving c into a strange shape different from 4N * 1, solver.error")

