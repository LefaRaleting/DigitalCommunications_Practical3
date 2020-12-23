# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 13:09:12 2020

@author: user
"""
import numpy as np
from itertools import combinations
from operator import add
import math
import matplotlib.pyplot as plt

#-----------------------------------------------------------------------------
#                           Varivables 
#-----------------------------------------------------------------------------
size= 100


G_1= [[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0],
     [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1],
     [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0],
     [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1],
     [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0],
     [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1],
     [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0],
     [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1],
     [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1],
     [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1]]

#G_1 = np.array(G_1).reshape((-1, 14)).tolist()
#s_1= BPSK(bitgenerator(size))
#s_1 = [1, 0, 0, 1, 0, 1, 1, 1]

#-----------------------------------------------------------------------------
#                           Bit generator 
#-----------------------------------------------------------------------------

def bitgenerator(size):
    
    return  [np.random.randint(2) for i in range(size)]

#-----------------------------------------------------------------------------
#                          BPSK generator 
#-----------------------------------------------------------------------------
def BPSK(bits):
    bpsk = []
    for k in bits:
        if k == 1:
            bpsk.append(1)
        else:
            bpsk.append(-1)
    return bpsk

#-----------------------------------------------------------------------------
#                           Linear block code
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#                          Parity check matrix (H)
#-----------------------------------------------------------------------------

def  parityCheck( G):
    G=np.array(G)
    p= G[0:G.shape[0],G.shape[0]:G.shape[1]] 
    
    H= np.concatenate((p.transpose(),np.identity(G.shape[1]-G.shape[0])),1)
    
    return H

#-----------------------------------------------------------------------------
#                          Encoding
#-----------------------------------------------------------------------------

def codeword(sent,G):#sent is a list , G is a list, n is noise
    G=np.array(G)
    s=np.array(sent).reshape((-1,G.shape[0]))
    codewords= s.dot(G) %2
    return codewords.flatten()


#-----------------------------------------------------------------------------
#                          Decoding
#-----------------------------------------------------------------------------
  
def  decoder(Recievedbits,G): #G is a list and Rec is a list too
    G=np.array(G)
    c=np.array(Recievedbits).reshape((-1,G.shape[1]))
    H=parityCheck(G)
    z=c.dot(H.T)%2
    
    for k in range(0,len(c)):
        if(sum(z[k]>=1)):
            index=[]
            for i in H.T: # Search for the syndrom Vec in H.T matrix
                if( (i == z[k]).all()):
                    index.append(1)
                else:
                    index.append(0)
            #check if z was found in the H.T matrix        
            if(sum(index)>0): #found something change  cest
                c[k] ^= np.array(index) #xor
                
            else: #nothing was found so find combinations
                for com in   combinations(range(len(H.T)),2):
                    OR = H.T[com[0]].astype(int)^H.T[com[1]].astype(int)
                    if(OR==z[k]).all():
                        c[k][com[0]] ^= 1
                        c[k][com[1]] ^= 1
                        break
    return c[:,:G.shape[0]].flatten()
#------------------------------------------------------------------------------
def Add_noise(transmitted,RC, SNR):
    M=2
    Gnoise= np.random.normal(0,size=len(transmitted))
    gama = 1 / np.sqrt(math.pow(10, (SNR / 10)) * 2 * RC)
    # print(gama)
    new = [i * gama for i in Gnoise]
    R = list(map(add, transmitted, new))
    return R

#k=[1,-1,1,1]


# ____________________________________________________________________________
#                           Detection
# ____________________________________________________________________________

def BPSKDetection(comp):
    points = [-1, 1] #sybols
    Bpoints = [[0], [1]]
    recieved = -1
    minDistance = 99
    decoded = []
    Bdecoded = []
    for y in comp:
        if(y>0):
            decoded.append(1)
            Bdecoded.append(1)
        else:  
            decoded.append(-1)
            Bdecoded.append(0)
            
    return decoded, Bdecoded  # recieved


#__________________________________________________________________________
#                      Bit Error calculation
#______________________________________________________________________________

def bit_errors(sent, recieved):
    error = 0
    for k in range(0,len(recieved)):
        if sent[k] != recieved[k]:
            error += 1
    BER = error / len(recieved)*100

    return BER



s_1= bitgenerator(100)
def tester():
    
    xValues = np.linspace(-4, 15, 100)
    yvalues=[]
    yvalues1=[]
    
    for x in xValues:
        ber=0
        ber1=0
        for i in range(0,10):
            s_1= bitgenerator(100)
            symbol,bits=BPSKDetection(Add_noise(BPSK(codeword(s_1,G_1).tolist()),0.5,x))
            symbol1,bits1=BPSKDetection(Add_noise(BPSK(s_1),1,x))
        #bits = [item for sublist in bits for item in sublist]
        #print(Add_noise(codeword(s_1,G_1).tolist(),0.5,1000))
        #print(bits)
        #c_est = [1,0,0,0,0,1,1,1,1,0,1,0,1,0] 
            ber+= bit_errors(s_1,decoder(bits,G_1).tolist())
            ber1+= bit_errors(s_1,bits1)
            
        ber=ber/10
        ber1=ber1/10
        yvalues.append(ber)
        yvalues1.append(ber1)
    plt.semilogy(xValues,yvalues, label="Linear Block Code")
    plt.ylabel('BER')
    plt.xlabel('SNR (dB)')
    
    plt.semilogy(xValues,yvalues1, label="Optimal detection uncoded")
    
    
    
    
    #At the end
    plt.title(" BER vs SNR")
    plt.legend()

tester() 