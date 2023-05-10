#Emma Klemets, Jan 2023
#Functions to calculate the adiabatic parameter

import numpy as np
import sympy as sp
import scipy.constants as const
from scipy import ndimage
from numpy import linalg as LA

#importing physical constants that are useful from scipy.constants
m_p = const.physical_constants['proton mass energy equivalent in MeV'][0]
m_n = const.physical_constants['neutron mass energy equivalent in MeV'][0]
c = const.c #m/s
hbar = const.hbar

mu_n = const.physical_constants['neutron mag. mom.'][0] #J T^-1
g_n = const.physical_constants['neutron mag. mom. to nuclear magneton ratio'][0]
#correction for frequency - should still check this with someone else
gamma_n = const.physical_constants['neutron gyromag. ratio'][0]/(2*np.pi) #s^-1 T^-1, negative


#3.11 in Pierre
def K_equ_dbdt(B, dbdt):
    return gamma_n*B**2 / dbdt

#equation 2 in the CDR - Sect 4.1 for the value of k
#3.12 in Pierre
def K_equ_dbdx(B, db_perpdx, v_n):
    return gamma_n*B**2 / (v_n*db_perpdx)

#3.14 in Pierre
def K_equ3_14(v_vec, B_1, B_2, r_1, r_2, k_inf_set=False):
    '''
    I trust this one a lot more
    equation 3.14 from https://tel.archives-ouvertes.fr/tel-00726870 by E. Pierre

    but this one really returns a k value for the range between two points, not actually
    a point itself.

    Problems: 
    - not totally sure what the just B should be
    - how small a step this is valid for?
    '''
    #what to use for this value of the field?
    B = LA.norm(B_1)
    
    v_n = LA.norm(v_vec)
    B_1_norm = LA.norm(B_1)
    B_2_norm = LA.norm(B_2)
    
    if k_inf_set:
        k_inf=np.inf
    else:
        k_inf=9e9        
    
    #the angle that the field changes by
    if B_1_norm == 0 or B_2_norm == 0: # if either 0 feild == bad, k=0
        print("Field is 0")
        k = 0
    elif B_1_norm*B_2_norm == 0: 
        '''
        This and the next case catch:
        if the change is ~0, that's fine, as it means the field is definitely
        changing slow enough for the spin to follow. So this is essentially 
        the k=inf case, but for graphing I'll just set it to a large value
        '''
        # print(B_1_norm, B_2_norm)
        k = k_inf
        
    elif B_1@B_2/(B_1_norm*B_2_norm) >= 1.0: 
        # print(B_1_norm, B_2_norm)
        k = k_inf
    else:
        theta = np.arccos(B_1@B_2/(B_1_norm*B_2_norm))
        # print(f"theta:{theta}, arg:{B_1@B_2/(B_1_norm*B_2_norm)}")

        #the distance for that change
        delta_distance = LA.norm(r_2 - r_1)

        # print(v_n, theta)
        #k calculation
        k = gamma_n*B*delta_distance / (v_n * theta)
    return k


#calculates the k value for a 90 degrees turn with a uniform magnitude magnetic field (B)
#for a given neutron speed (vn) and distance of the turn (dx) 
def pseudoScalar_K_equ3_14_turn(vn, B, dx, k_inf_set=False):
    
    v_vec = np.array([vn, 0, 0])
    B_1 = np.array([B, 0, 0])
    B_2 = np.array([0, B, 0])
    r_1 = np.array([0, 0, 0]) 
    r_2 = np.array([dx, 0, 0])

    return K_equ3_14(v_vec, B_1, B_2, r_1, r_2, k_inf_set=k_inf_set)


#calculates the k value a field that ic changing, but assuming the neutron
#is moving stright along one axis 
#for a given neutron speed (vn)
def pseudoScalar_K_equ3_14(vn, B, x, k_inf_set=False):
    K_arr = []
    
    for i in range(len(B)-1):
    
        v_vec = np.array([vn, 0, 0])
        B_1 = np.array([B[i], 0, 0])
        B_2 = np.array([0, x[i+1], 0])

        r_1 = np.array([x[i], 0, 0]) 
        r_2 = np.array([x[i+1], 0, 0])
        
        
        K_arr.append(K_equ3_14(v_vec, B_1, B_2, r_1, r_2, k_inf_set=k_inf_set))

    return K_arr


#from just the B field
def K_equ14(v_n, B, dB_vec_dy, k_inf_set=False):

    B_norm = LA.norm(B)

    denom  = LA.norm(np.cross(dB_vec_dy, B))

    k = gamma_n*B_norm**3 / (v_n * denom)

    return k
