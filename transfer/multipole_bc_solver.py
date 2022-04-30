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

from calendar import c
from typing import Any
import numpy as np
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
def get_modal_displacement(file)->Any:
    'This function should read in the modal displacement'
    'of M Boundary Conditions[BCs]'
    "这是一个用于读入“模态边界振动情况”的函数"
    ## file_in operation
    ## read operation
    modal_displacement = ... file 
    return modal_displacement 
    ## this modal_displacement should be a M * 1*3 vector?
    ## 返回的应该是有M个1*3向量的数组
    
def get_boundary_normal_vector(file)->Any:
    'This function should read in the normal vector of the boundary of the body'
    "这是一个读入物体轮廓法向的函数"
    normal_vector = ... 
    return normal_vector
    ## the normal_vecotr should be a M * 1 * 3 vector
    ## 返回的应该是有M个1*3向量的数组
    
def construct_rhs(normal_vector, modal_vector)-> Any:
    
    for i in range(0, M)

    return 
    