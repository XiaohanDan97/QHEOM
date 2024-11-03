# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 15:49:58 2024

@author: Xiaohan Dan
"""

import numpy as np
import scipy.linalg as LA
from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit, transpile
from qiskit.quantum_info.operators import Operator
import walsh_gray_optimization as wo

#generate the quantum gate matrices for SVD-dilation
def gene_umat_SVD(Gmat,Isscale=True):
  #Dilation process: Using SVD.
  U1,S1,V1 = LA.svd(Gmat)
  #large normfac will lead to less good result
  if(Isscale):
    normfac = S1[0] #*1.1
    S1 = S1/normfac
  else:
    normfac = 1.0
  Nvec = Gmat.shape[0]

  #Numerically, norm of Gprop should <=1 (<1 is best).
  Mzero = np.zeros((Nvec,Nvec),dtype=np.complex128)
  fk=np.zeros(2*Nvec,dtype=np.float_)

  Sig_p = np.zeros(Nvec,dtype=np.complex128)
  Sig_m = np.zeros(Nvec,dtype=np.complex128)
  for i in range(len(S1)):
    Sig_p[i] = S1[i]+1j*np.sqrt((1-S1[i]**2))
    Sig_m[i] = S1[i]-1j*np.sqrt((1-S1[i]**2))

    #here U = e^{i f_k}
    fk[i] = (-1j*np.log(Sig_p[i])).real
    fk[Nvec+i] = (-1j*np.log(Sig_m[i])).real

  SG0 = np.block([[np.diag(Sig_p),Mzero],\
        [Mzero,np.diag(Sig_m)]])
  
  return normfac,U1,V1,SG0,fk

#for re-calculating normfac_list and time_list
def cal_tnorm_list(timearr,G,Nlen,Isscale=True):
  normfac_list = []
  time_q = np.zeros(Nlen)
  nlarge = int(len(timearr)/(Nlen-1))
  for i0 in range(Nlen):
    istep = i0*nlarge
    time_q[i0] = timearr[istep]
    normfac,U_matu,U_matv,S_mat0,fk = gene_umat_SVD(G[istep],Isscale)
    normfac_list.append(normfac)
  return time_q, normfac_list

#Construct the SVD+Walsh circuit
def SVD_walsh_cirq(Nqb, Gmat, ini_vec):

  qc = QuantumCircuit(Nqb,Nqb)

  normfac,U_matu,U_matv,S_mat0,fk = gene_umat_SVD(Gmat)

  #the walsh coeff
  arr_a = wo.walsh_coef(fk,Nqb)
  qc.append(Operator(U_matv),range(0,Nqb-1))
  qc.h(Nqb-1)

  Ulist_diag0 = wo.cirq_list_walsh(arr_a,Nqb,1E-5)
  Ulist_diag = wo.optimize(Ulist_diag0)
  qc_diag = wo.cirq_from_U(Ulist_diag,Nqb)

  qc.append(qc_diag.to_gate(),range(Nqb))
  qc.append(Operator(U_matu),range(0,Nqb-1))

  qc.h(Nqb-1)

  return qc,normfac

#Construct the SVD-dilation circuit (without Walsh Operator)
def dilation_SVD_cirq(Nqb,Gt,ini_vec):

  normfac,U_matu,U_matv,S_mat,fk = gene_umat_SVD(Gt)

  #construct a circuit, (num of qubits, num of classical bits)
  circuit=QuantumCircuit(Nqb,Nqb)
  circuit.initialize(ini_vec,range(0,Nqb))
  
  circuit.append(Operator(U_matv),range(0,Nqb-1))
  circuit.h(Nqb-1)
  circuit.append(Operator(S_mat),range(Nqb))
  circuit.append(Operator(U_matu),range(0,Nqb-1))
  circuit.h(Nqb-1)

  return circuit,normfac

#Construct the Traditional Sz-Nagy circuit (without Walsh Operator)
def Sz_Nagi_cirq(Nqb,array,ini_vec):

  #construct a circuit, (num of qubits, num of classical bits)
  circuit=QuantumCircuit(Nqb,Nqb)
  circuit.initialize(ini_vec,range(0,Nqb))

  # Normalization factor, 1.1 times martix's norm to ensure contraction
  norm = LA.norm(array,2)*1.1
  array_new = array/norm

  ident = np.eye(array.shape[0])

  # Calculate the conjugate transpose of the G propagator
  fcon = (array_new.conjugate()).T

  # Calculate the defect matrix for dilation
  fdef = LA.sqrtm(ident - np.dot(fcon, array_new))

  # Calculate the defect matrix for the conjugate of the G propagator
  fcondef = LA.sqrtm(ident - np.dot(array_new, fcon))

  # Dilate the G propagator to create a unitary operator
  array_dilated = np.block([[array_new, fcondef], [fdef, -fcon]])

  circuit.append(Operator(array_dilated),range(Nqb))

  return circuit, norm