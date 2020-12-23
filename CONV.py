# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 17:56:06 2020

@author: user
"""

import numpy as np
from copy import deepcopy
from collections import deque
import matplotlib.pyplot as plt

from LBC import bitgenerator,BPSK,Add_noise,BPSKDetection,bit_errors,decoder,decoder

K=3

class MNode:
    def __init__(self, status):
        self.status = status
        self.heuristic = None
        self.Predecessor = None


class Viterbi_MLSE:
    def __init__(self, received, G):
        self.K = len(G)
        self.received = np.array(received).reshape((-1, self.K))
        self.BlockN = len(self.received)
        self.Root = None
        
    def Surving_path(self, lists, new_item):
        for item in lists:
            if new_item.status == item.status:
                if new_item.heuristic > item.heuristic:
                    return False
                else:
                    lists.remove(item)
                    lists.append(new_item)
                    return True
                    
        lists.append(new_item)
        return True
    
    
    def Gen( self,val):
        val.reverse()
        c1= val[0]
        c2= val[0]^val[2]
        c3= val[0]^val[1]
        
        return [c1,c2,c3]

    def ShiftRegister(self, x_state, m):
        status = deepcopy(x_state)
        status.pop()
        status = list(reversed(status))
        status.append(m)
        return list(reversed(status))


    

    def MLSE_Decode(self):
        bits = [1, 0]
        self.Root =  MNode([0,0] )
        self.Root.heuristic=0.0

        New_list = []# deque()
        New_list.append(self.Root)
        tail = None

        
        for t in range(self.BlockN):
            Old_list = deepcopy(New_list)
            New_list.clear()
            
            if t <= self.BlockN - (self.K - 1):
                while len(Old_list) > 0:
                    T_root = Old_list.pop()
                    for bit in bits:
                        T_node = MNode(self.ShiftRegister(T_root.status, bit))
                        T_node.Predecessor = T_root
                        heuristic = self.calculate_cost(self.received[t],
                                                   self.Gen(T_node.status + [T_root.status[1]]))
                        
                        T_node.heuristic = heuristic + T_root.heuristic

                        self.Surving_path(New_list, T_node)
            else:
                while len(Old_list) > 0:
                    T_root = Old_list.pop()
                    # shift in 0
                    T_node = MNode(self.ShiftRegister(T_root.status, 0))
                    T_node.Predecessor = T_root
                    heuristic = self.calculate_cost(self.received[t],
                                               self.Gen(T_node.status + [T_root.status[1]]))
                   
                    T_node.heuristic = heuristic + T_root.heuristic

                    if self.Surving_path(New_list, T_node):
                        tail = T_root
            
           
           

        temp_tail = tail
        arr1 = []
       
        
        while temp_tail is not None:
            arr1.append(temp_tail.status[0])
            temp_tail = temp_tail.Predecessor
        send = list(reversed(arr1))
        
        return [0] +send + [0]

    def calculate_cost(self, conv, conv_est):
        conv_est = 2 * np.array(conv_est) - 1
        #conv = 2* np.array(conv)-1
        delta = list(map(abs,conv - conv_est))
        delta = np.power(np.array(delta),2)
        return delta.sum()

    






#------------------------------------------------------------------------------


def Gen( val):
    val.reverse()
    c1= val[0]
    c2= val[0]^val[2]
    c3= val[0]^val[1]
    
    return [c1,c2,c3]

def encoder(sentbits): #this function is used for encoding 
    codeword=[]
    for k in range(len(sentbits)-2):
        val= sentbits[k:k+3] #reads 3 bits at a time
        codeword+=Gen(val)
        
    return codeword


s_1= bitgenerator(100)
def tester():
    
    xValues = np.linspace(-4, 15, 30)
    yvalues=[]
    yvalues1=[]
    
    for x in xValues:
        ber=0
        ber1=0
        for i in range(0,10):
            s_1= bitgenerator(100)
            symbol,bits=BPSKDetection(Add_noise(BPSK(encoder(s_1)),0.5,x))
            symbol1,bits1=BPSKDetection(Add_noise(BPSK(s_1),1,x))
            Q= Viterbi_MLSE(symbol,[4,5,6])
            
            ber+= bit_errors(s_1,Q.MLSE_Decode())
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