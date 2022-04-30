
from tabnanny import check
from typing import Any
import numpy as np
import multipole_bc_solver

# macro
PI = np.pi
OMEGA = 2763    # a specified frequence, in Hz, will be a constant if the modal is specified
c = 340         # in m/s velocity of sound in the air.
k = OMEGA/c

class Multipole:
    coeffiency = [0, 0, 0, 0]
    loc = [0, 0, 0]
    def __init__(self, coeffiency, loc) -> None:
        self.coeffiency = coeffiency
        self.loc = loc
        pass
    # def __getattribute__(self, __name: str) -> Any:
    #     pass
    
def S_0_0(k, listener_pos, multipole_loc):
    "The S_0_0 function takes into multipole_loc[1 * 3] and listener_pos [1 * 3] as input and an pressure field [1 * 1] as output"
    r = np.sqrt(np.sum((np.array(listener_pos) - np.array(multipole_loc)) ** 2))
    res = 1/2 * np.sqrt(1/PI) * (1.0j/(k * r)) * np.exp(-1.0j*k*r)
    return res

def S_1_minus1(k, listener_pos, multipole_loc):
    "The S_1_-1 function takes into multipole_loc[1* 3] and listener_pos [1 * 3] as input and an pressure field [1 * 1] as output"
    r = np.sqrt(np.sum((np.array(listener_pos) - np.array(multipole_loc)) ** 2))
    diff = np.array(listener_pos) - np.array(multipole_loc)
    diff_x = diff[0]
    diff_y = diff[1]
    res = -1/2* np.sqrt(3/(2 * PI)) * ((diff_x - 1.0j* diff_y) / r ) * ((k * r - 1.0j) / ((k * r) ** 2)) * np.exp(-1.0j * k * r)
    return res

def S_1_0(k, listener_pos, multipole_loc):
    "The S_1_-1 function takes into multipole_loc[1* 3] and listener_pos [1 * 3] as input and an pressure field [1 * 1] as output"
    r = np.sqrt(np.sum((np.array(listener_pos) - np.array(multipole_loc)) ** 2))
    diff = np.array(listener_pos) - np.array(multipole_loc)
    diff_z = diff[2]
    res = -1/2* np.sqrt(3/(2 * PI)) * (diff_z / r ) * ((k * r - 1.0j) / ((k * r) ** 2)) * np.exp(-1.0j * k * r)
    return res



def S_1_1(k, listener_pos, multipole_loc):
    "The S_1_-1 function takes into multipole_loc[1* 3] and listener_pos [1 * 3] as input and an pressure field [1 * 1] as output"
    r = np.sqrt(np.sum((np.array(listener_pos) - np.array(multipole_loc)) ** 2))
    diff = np.array(listener_pos) - np.array(multipole_loc)
    diff_x = diff[0]
    diff_y = diff[1]
    res = -1/2* np.sqrt(3/(2 * PI)) * ((diff_x + 1.0j* diff_y) / r ) * ((k * r - 1.0j) / ((k * r) ** 2)) * np.exp(-1.0j * k * r)
    return res

def  get_multipole_coefficience_series(coefficiency_file)->Any:
    'This function should read in the coefficience of N multipole'
    "这是一个读入4N个多级源系数的函数"
    ...
    coefficiency_file
    ...
    coefficiency_series = ...
    return coefficiency_series
    ## the multipole_coefficiency_series should be a 4N*1 vector
    ## 直接返回4N*1的向量就可以了
    
def  get_listening_position(listen_pos_file)->Any:
    'This function should read in the coefficience of N multipole'
    "这是一个读入4N个多级源系数的函数"
    ...
    listen_pos_file
    ...
    listen_pos_series = ...
    return listen_pos_series
    ## the multipole_location_series should be a K * 1 * 3 vector
    ## 返回的应该是有K个1*3向量的数组
    
    
def check_in(coefficiency_file,multipole_location_file):
    "This function takes read in c[4N * 1] and multipole_loc_sr[N * 1 * 3]"
    multipole_loc_sr = multipole_bc_solver.get_multipole_location_series(multipole_location_file=multipole_location_file)
    c_sr = get_multipole_coefficience_series(coefficiency_file=coefficiency_file)
    if len(c_sr) / 4 == len(multipole_loc_sr):
        N = len(multipole_loc_sr)
        return N, c_sr, multipole_loc_sr
    else:
        print("while preparing rendering sound, the multipole coefficiency length does not match the location length, check_in error")

def pressure_computing(pos, k, c_sr, multipole_loc_sr):                                                  # we still need file here?
    "This function takes into c[4N * 1] and multipole_loc_sr[N * 1*3] the listening position[1*3] and output the pressure in this position[1*1]"
    '读入4N*1的系数和N*1*3的多级源位置，再读入波数和听者位置，给出听者位置的声压'
    if len(c_sr) / 4 == len(multipole_loc_sr):
        N = len(multipole_loc_sr)
    else:
        return 
    pressure = 0
    for i in range(0, N):
        if i % 4 == 0:
            pressure += c_sr[i] * S_0_0(k=k, listener_pos=pos, multipole_loc=multipole_loc_sr[i])
        if i % 4 == 1:
            pressure += c_sr[i] * S_1_minus1(k=k, listener_pos=pos, multipole_loc=multipole_loc_sr[i])
        if i % 4 == 2:
            pressure += c_sr[i] * S_1_0(k=k, listener_pos=pos, multipole_loc=multipole_loc_sr[i])
        if i % 4 == 3:
            pressure += c_sr[i] * S_1_1(k=k, listener_pos=pos, multipole_loc=multipole_loc_sr[i])
    return pressure

def multipoint_pressure_computing(listen_pos_series, k, c_sr, multipole_loc_sr):
    "This function into c[4N * 1] and multipole_loc_sr[N * 1*3] the listening position[K * 1*3] and output the pressure in this position[K * 1*1]"
    K = len(listen_pos_series)
    pressure_series = np.zeros((K, 1))
    for i in range(0, K):
        pressure_series[i] = pressure_computing(pos=listen_pos_series[i], k=k,c_sr=c_sr, multipole_loc_sr=multipole_loc_sr)
    return pressure_series

def sound_render(listen_pos_file, k, coefficiency_file, multipole_location_file):
    N, c_sr, multipole_location_sr = check_in(coefficiency_file=coefficiency_file, multipole_location_file=multipole_location_file)
    listen_pos_sr = get_listening_position(listen_pos_file=listen_pos_file)
    pressure_sr = multipoint_pressure_computing(listen_pos_series=listen_pos_sr, k=k, c_sr=c_sr, multipole_loc_sr=multipole_location_sr)
    
    fileout 
    pressure_sr
    return pressure_sr