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

from audioop import mul
from calendar import c
from cmath import pi
from typing import Any
import numpy as np
from pandas import array
from a_scratch import OMEGA
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
def get_modal_displacement_series(file)->Any:
    'This function should read in the modal displacement'
    'of M Boundary Conditions[BCs]'
    "这是一个用于读入“模态边界振动情况”的函数"
    ... ## file_in operation
    file
    ... ## read operation
    modal_displacement_series = ... 
    return modal_displacement_series 
    ## this modal_displacement should be a M * 1*3 vector?
    ## 返回的应该是有M个1*3向量的数组
    
def get_boundary_normal_vector_series(file)->Any:
    'This function should read in the normal vector of the boundary of the body'
    "这是一个读入物体轮廓法向的函数"
    ...
    file
    ...
    normal_vector_series = ... 
    return normal_vector_series
    ## the normal_vecotr should be a M * 1 * 3 vector
    ## 返回的应该是有M个1*3向量的数组

def  get_multipole_location(file)->Any:
    'This function should read in the location of M multipole'
    
    
def construct_rhs(normal_vector, modal_displacement)-> Any:
    "construct right hand side"
    "This function takes in normal_vector and modal_displacement"
    "gives out the boundary condition"
    '构造方程右边，读入M * 1 * 3的“法向向量”和M * 1 *3的“模态振动？？”，输出M*1的边界条件bc'
    if(len(normal_vector) == len(modal_displacement)):
        M = normal_vector
    else:
        print("the size of normal_vecor and the modal_vector do not match, construct_rhs error")
    a = np.array(normal_vector)
    b = np.array(modal_displacement)
    c = (a * b).sum(axis = 1)
    bc = OMEGA ** 2 * RHO * c           ## bc means boundary_condition
    if(np.shape(bc)==(M,)):
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
    res = -1/2 * np.sqrt(1/pi) * (1/(k * r**3)) * (0.0 + 1.0j) * np.exp(0.0 - 1.0j * k * r) * (1.0 + 1.0j * k * r) * delta
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
    res = -1/2 * np.sqrt(3 / (2 * pi)) * np.exp(0.0 - 1.0j*k*r) * (1/(k**2*r**5)) * ((k * r - 1.0j)*delta1 + (-2*k*r + 1.0j * (2 - k**2*r**2))*(diff_x-1.0j*diff_y)*delta0)   
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
    res = -1/2 * np.sqrt(3 / (2 * pi)) * np.exp(0.0 - 1.0j*k*r) * (1/(k**2*r**5)) * ((k * r - 1.0j)*delta1 + (-2*k*r + 1.0j * (2 - k**2*r**2))*(diff_z)*delta0)   
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
    res = -1/2 * np.sqrt(3 / (2 * pi)) * np.exp(0.0 - 1.0j*k*r) * (1/(k**2*r**5)) * ((k * r - 1.0j)*delta1 + (-2*k*r + 1.0j * (2 - k**2*r**2))*(diff_x+1.0j*diff_y)*delta0)   
    return res

def construct_coefficiency_matrix(sampling_loc, multipole_loc, normal_vector)->Any:
    
    return