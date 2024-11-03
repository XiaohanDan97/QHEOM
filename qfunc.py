import numpy as np
from qiskit import QuantumRegister, QuantumCircuit

#read and write
X = np.array([[0,1],[1,0]])
Y = np.array([[0,-1j],[1j,0]])
Z = np.array([[1,0],[0,-1]])
I = np.array([[1,0],[0,1]])

def read_superoper_array(TIME_STEPS,Para):
    S_array_read = np.zeros((TIME_STEPS, Para.Nmat, Para.Nmat), dtype=np.complex_)
    time = np.zeros((TIME_STEPS), dtype=np.float_)
    for i in range(Para.Nmat):
        for j in range(Para.Nmat):
            time_read, Sreal, Simag = np.hsplit(np.loadtxt("G_"+str(i+1)+str(j+1)+ ".dat"), 3)

            # stores values in variables and combines real and imag parts of super-operator
            for k in range(TIME_STEPS):
                time[k] = time_read[k][0]
                S_array_read[k][i][j] = Sreal[k][0] + 1.j * Simag[k][0]
    return time,S_array_read

#expand the propagator
def expand(S_arr,dim1):
    TIME_STEPS,dim0,dim02 = S_arr.shape
    if(dim0!=dim02): print('err in expand 1')
    S_arr_exp = np.zeros((TIME_STEPS, dim1, dim1),dtype=np.complex_)
    S_arr_exp[:TIME_STEPS,:dim0,:dim02] = S_arr
    for i in range(dim0,dim1):
        for istep in range(TIME_STEPS):
            S_arr_exp[istep,i,i]=1.0
    return S_arr_exp

#project the propagator
#proj_state_list: the state needed
#0: population D, 1: Coherence, 2 Coherence, 3: population A
def proj(Gt,proj_state_list):
  nlen = len(proj_state_list)
  time_step = Gt.shape[0]
  Gt_new = np.zeros((time_step, nlen, nlen), dtype=np.complex_)
  for i in range(time_step):
    for j in range(nlen):
      for k in range(nlen):
        Gt_new[i,j,k] = Gt[i,proj_state_list[j],proj_state_list[k]]
  return Gt_new
