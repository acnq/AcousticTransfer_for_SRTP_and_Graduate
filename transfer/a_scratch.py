
from asyncio import streams
from time import time
from tkinter import N
from webbrowser import get
import scipy.linalg
import numpy as np

## macro:
RHO = 1.2041    # density of air in kg/m^3, for compute the boundary condition
OMEGA = 2763       # specified frequence in Hz, will be a constant if modal specified 
c = 340         # in m/s velocity of sound in the air, 


# function1:
# construct_rhs: 
#       construct the Right Hand Side of linear_sys, we need 
# u_sam_l:
#       the position of l-th particle after vibration in the discreted deformed body
# u_sam_0:
#       the equilibrium position of lth particle in the  discreted deformed body
# n_x:
#       the normal vector of the surface surronding the l-th particle 
# omega:
#       the frequency of the specified modal.
def construct_rhs(u_sam_l, u_sam_0, n_x, omega=OMEGA):
    d_sam = u_sam_l - u_sam_0
    a = (d_sam * n_x).sum(axis=1) 
    # python中*表示各元素相乘，利用.sum(axis=1)表示横向求和
    rhs = RHO * a * omega ** 2
    return rhs

# 以下代码需要重新测试！！！（TODO）
# functions2:
# construct_matrix:
#       construct the coefficiency matrix of the linear_sys, we need:
# n_x:
#       the vector of the surface of the body
# N:
#       the number of the multipole you want to use
# u_x_l:
#       the position of sampling points
# u_mul:
#       the position of multipole expansion position

def construct_matrix(N, n_x, u_mul, u_sam,):
    '相关下标注意重新再修正一下，可能依旧存在问题'
    A = np.zeros((N, N))
    R = np.zeros((N, N, 3))
    # print("R:")
    for l in range(0, N):
        for m in range(0, N):
            for k in range(0, 3):
                R[m, l, k] = u_mul[m, k] - u_sam[l, k]
    # print(R)
    
    # print("dist")
    dist = np.zeros((N, N))
    for l in range(0, N):
        for m in range(0, N):
            dist[l, m]=np.sqrt((u_mul[l,0]- u_sam[m,0]) ** 2 + (u_mul[l,1] - u_sam[m,1]) ** 2 + (u_mul[l, 2] - u_sam[m, 2]) ** 2)
    # print(dist)

    K = OMEGA / c
    A1 = np.zeros((N, N), dtype=complex)
    A2 = np.zeros((N, N), dtype=complex)
    A = np.zeros((N, N), dtype=complex)
    for l in range(0, N): # l 指标sampling points
        for m in range(0, N): ## m指标multipole points,代替推导中的j，似乎python中已经用j代表虚数单位了
            A1[l, m] = np.exp(k * dist[l,m] * (.0 - 1.0j)) * (1./(dist[l, m] ** 2) - (0 + 1.0j)/(K ** 2 * dist[l,m] ** 2))
            A2[l, m] = (R[m, l, 0] * n_x[l, 0] + R[m, l, 1] * n_x[l, 1] + R[m, l, 2] * n_x[l, 2]) / dist[m, l] #注意n的指标是l,和sampling point 的指标一样
            A[l, m] = A1[l, m] * A2[l, m]
    # print("A1: ")
    # print(A1)
    # print("A2:")
    # print(A2)
    # print("A:")
    # print(A)

    return A

u_x_0 = np.array([[1.324, 2.754, 3.758], [5.589, 7.645, 9.137], [7.324, 8.234, 9.532]])               # 用户提供的采样点初始位置
u_x_l = np.array([[1,2,3], [5,7,9], [7,8,9]])               # 用户提供的采样点位移后的位置
## 我们通过以上两种信息去确定用户提供的采样点的“位移”如果用户能够直接提供采样点的位移，则可以直接提供
## 我这边修改一下函数的接口就可以完成了
n_x = np.array([[.1, .2, .3], [.4, .5, .6], [.7,.8,.9]])    # 用户提供的物体轮廓的法向方向

u_sam=np.array([[1,2,3], [5,7,9], [7,8,9]])                   ## 这个应该是用户提供的采样点的初始位置
u_mul=np.array([[4,5,9], [3,7,6], [5, 8, 1]])                                    ## 这个应该是用户提供的采样点位置


## 总之，目前这套程序至少需要用户提供：1.采样点位置， 2. 采样点位移，3. 采样点所在的物体的边界的法向，4. 多极源展开点的位置
rhs = construct_rhs(u_x_l, u_x_0, n_x)
matrix = construct_matrix(3, n_x, u_mul=u_mul, u_sam=u_sam)

print("rhs:", rhs)
print("matrix:" ,matrix)

x = np.linalg.solve(matrix, rhs) ## 今天单纯使用solver, 之后可能要调一下np.lisqrt
print("solution1:", x)


p, res, rnk, s = scipy.linalg.lstsq(matrix, rhs) ## 似乎直接写成scipy.linalg.lstsq会出问题
print("solution2:", p) 

# to create a sound, 
# one should specified the position of listener(x)
# the position of multipoles (u_mul)
# and the "now determined coefficients p"

x = np.array([1000, 1100, 1200])
def sound_pressure(N, p, x, u_mul, omega=OMEGA, c=c):
    R = np.zeros(N)
    for l in range(0, N):
        R[l] = np.linalg.norm(x - u_mul[l])

    
    k = OMEGA / c
    K = np.zeros(N)
    for l in range(0, N):
        K[l] = (.0 + 1.0j) / (k * R[l]) * np.exp((.0 - 1.0j) * (k * R[l]))
    sound_pressure = np.dot(p, K)
    return sound_pressure

sound_pressure = sound_pressure(3, p, x, u_mul=u_mul)
print("sound pressure:", sound_pressure)

import pyaudio
duration = 1          # time period for sound play, in seconds 
fs = 44100              # sampling rate for a sound in Hz                
def get_sound(sound_pressure, fs, omega=OMEGA):
    ## sound = sound_pressure(x) * q(t)
    time = np.arange(omega * duration) / omega
    samples = (np.sin(omega * time) * sound_pressure ).astype(np.float32)
    
    p = pyaudio.PyAudio()
    stream = p.open(format=pyaudio.paFloat32,
                    channels=1,
                    rate=omega,
                    output=True)
    print(samples)
    # play
    stream.write(samples.tobytes())
    stream.stop_stream()
    stream.close()
    p.terminate()
    
    # save to wav
    from scipy.io.wavfile import write
    write('./transfer/result/result.wav', fs, samples)
    return samples
    
# def save_sound(samples):
#     # save to wav 2
#     import wave
#     WAVE_OUTPUT_FILENAME="output.wav"
#     frames=[]
#     wf = wave.open(WAVE_OUTPUT_FILENAME, 'wb')
#     wf.setnchannels(1)
#     wf.setsampwidth(p.get_sample_size(pyaudio.paFloat32))
#     wf.setframerate(fs)
#     wf.writeframes(b''.join(frames))
#     wf.close()


samples = get_sound(sound_pressure=sound_pressure, fs=fs)
# save_sound(samples=samples)