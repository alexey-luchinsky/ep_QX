import matplotlib.pyplot as plt
import numpy as np

def load_data(var,channel,dir='../build/'):
    return np.loadtxt(dir+'h'+var+'_'+channel+'.hst')

def load_exp(var,dir='../experiment/'):
    return np.loadtxt(dir+var+'.txt')

colList=['3S1_cs','3S1_co','1S0_co','3P0_co'];

# LDME
mc=1.5;
LDME={};
def set_LDME(r):
    LDME['3S1_cs']=0.270;
    LDME['3S1_co']=6.6e-3;
    LDME['1S0_co']=2.2e-2*3/(3*r + 1);
    LDME['3P0_co']=2.2e-2*3*mc**2*r/(3*r + 1);

def add(A,B):
    return np.array([ [A[i,0],A[i,1]+B[i,1]] for i in range(len(A))])
def mult(A,B):
    return np.array([ [A[i,0],A[i,1]*B[i,1]] for i in range(len(A))])
def div(A,B):
    return np.array([ [A[i,0],A[i,1]/B[i,1]] for i in range(len(A))])
def multC(c, A):
    return np.array( [ [A[i,0],c*A[i,1]] for i in range(len(A))])

data={};
def sum_th(var,r,dir='../build/'):
    set_LDME(r);
    all=multC(0,load_data(var,'3S1_cs',dir=dir))
    for ch in colList:
        data[ch]=multC(LDME[ch],load_data(var,ch,dir=dir))
        all = add(all,data[ch])
    return all
(rMin, rMax)=(0,1)

def get_xlabel(var):
    if var=='PT2':
        return r"$p_T^2,\.\mathrm{GeV}^2$"
    elif var=='Q2':
        return r"$Q^2,\.\mathrm{GeV}^2$"
    elif var=='W':
        return r"$W,\.\mathrm{GeV}$"
    else:
        return "???"

def get_ylabel(var):
    if var=='PT2':
        return r"$d\sigma/dp_T^2,\.\mathrm{pb}/\mathrm{GeV}^2$"
    elif var=='Q2':
        return r"$d\sigma/dQ^2,\.\mathrm{pb}/\mathrm{GeV}^2$"
    elif var=='W':
        return r"$d\sigma/dW,\.\mathrm{pb}/\mathrm{GeV}$"
    else:
        return "???"
