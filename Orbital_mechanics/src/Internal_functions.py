
import numpy as np
from Internal_constants import *

# --- FONCTIONS DE SIMULATION ---

def get_keplerian_acceleration(pos, M_central, G):
    """
    Calcul pur de l'accélération gravitationnelle vectorielle.
    Formule : a = - (G * M * pos) / |pos|^3
    """
    r_mag = np.linalg.norm(pos)
    # Protection contre la division par zéro
    if r_mag == 0: return np.zeros(3) 
    return - (G * M_central * pos) / r_mag**3

def compute_accel_nBody(R, Masses, nBody):
    # 1. Extraction et remise en forme (Phase Space)
    pos = R[:3*nBody].reshape((nBody, 3))
    vel = R[3*nBody:].reshape((nBody, 3))
    accel = np.zeros((nBody, 3))
    
    # 2. Somme des interactions par paires
    for i in range(nBody):
        for j in range(i + 1, nBody):
            # Vecteur de i vers j
            r_ij = pos[j] - pos[i]
            
            # Accélération de i subie par j (vers j)
            # On inverse le signe car get_keplerian renvoie un vecteur vers l'origine
            a_i_par_j = -get_keplerian_acceleration(r_ij, Masses[j], G_norm)
            
            # Accélération de j subie par i (vers i)
            a_j_par_i = -get_keplerian_acceleration(-r_ij, Masses[i], G_norm)
            
            accel[i] += a_i_par_j
            accel[j] += a_j_par_i
            
    # 3. Recombinaison du vecteur d'état dérivé
    return np.concatenate([vel.flatten(), accel.flatten()])


# --- PHYSIQUE : DÉRIVÉE POUR 2 CORPS (SOLEIL FIXE) ---
def compute_accel_2body(R_state, Masses):
    """
    R_state : [x, y, z, vx, vy, vz]
    Utilise la brique de base get_keplerian_acceleration
    """
    pos = R_state[:3]
    vel = R_state[3:]
    
    # Appel de la fonction atomique
    accel = get_keplerian_acceleration(pos, Masses[0], G_norm)
    
    return np.concatenate([vel, accel])

# --- MOTEURS RK4 ---

def rk4_step_nBody(R, Masses, nBody, h):
    """Pas RK4 pour le système complet (N corps)"""
    k1 = h * compute_accel_nBody(R, Masses, nBody)
    k2 = h * compute_accel_nBody(R + 0.5 * k1, Masses, nBody)
    k3 = h * compute_accel_nBody(R + 0.5 * k2, Masses, nBody)
    k4 = h * compute_accel_nBody(R + k3, Masses, nBody)
    
    return R + (k1 + 2*k2 + 2*k3 + k4) / 6

def rk4_step_2body(R, Masses, h):
    """Pas RK4 pour le système simplifié (Soleil fixe)"""
    k1 = h * compute_accel_2body(R, Masses)
    k2 = h * compute_accel_2body(R + 0.5 * k1, Masses)
    k3 = h * compute_accel_2body(R + 0.5 * k2, Masses)
    k4 = h * compute_accel_2body(R + k3, Masses)
    
    return R + (k1 + 2*k2 + 2*k3 + k4) / 6






# Protoypes pour dormand_prince, runge_kutta et euler:
# y: Vecteur d'état. (éléments Képlérien, x,y,z, )
# t: Epoque actuelle
# h: Pas d'intégration
# f: Fonction dans laquelle il y a les calculs pour intégrer les dérivées des variables d'état

"""
def dormand_prince(y, t, h, f):
    
    k1= h * f(t,y)
    k2= h * f(t+h/5, y + k1/5 )
    k3= h * f(t+3/10*h, y + 3/40*k1+ 9/40*k2)
    k4= h * f(t+4/5*h , y +44/45*k1 - 56/15*k2 + 32/9*k3)
    k5= h * f(t+8/9*h , y + 19372/6561*k1 - 25360/2187*k2 + 64448/6561*k3 - 212/729*k4)
    k6= h * f(t+h , y + 9017/3168*k1 - 355/33*k2 - 46732/5247*k3 + 49/176*k4 - 5103/18656*k5)
    k7= h * f(t+h , y + 35/384*k1 + 500/1113*k3 + 125/192*k4 - 2187/6784*k5 + 11/84*k6)
    
    y= y + 35/384*k1 + 500/1113*k3 + 125/192*k4 - 2187/6784*k5 + 11/84*k6
    return y

def runge_kutta(y, t, h, f):
    
    l1=h*f(t,y)
    l2=h*f(t+h/2 , y+ l1/2)
    l3=h*f(t+h/2, y + l2/2)
    l4=h*f(t+h,y+l3)
    
    y=y + l1/6 + 2/6*l2 +2/6*l3 + l4/6
    return y
   
def euler(y, t, h, f):
    y=y + h * f(t,y)
    return y


def dR(R,t):

    # Elements orbitaux

    a       = R[0]   
    e       = R[1]        
    i       = R[2]    
    OMEGA   = R[3]    
    w       = R[4]    
    M       = R[5]
    
    # expression de n
    
    n = 2 * np.pi
    
    # Derivées partielles du potentiel R qui s'arrête à J2
    
    dRda = -(9/2)*((n**2 * rt**2 * J2)/(a*(1-e**2)**(3/2))) * (1/3 - (np.sin(i)**2)/2)
    dRde = (9/2)*((n**2 * rt**2 * J2 * e)/((1-e**2)**(5/2))) * (1/3 - (np.sin(i)**2)/2)  
    dRdi= - (3/2)* ( n**2 * rt**2 * J2)/((1-e**2)**(3/2)) * np.cos(i)*np.sin(i)
    
    dRdw=0
    dRdOMEGA=0
    dRdM=0
    
    # système d'equations differentielles de Lagrange: da/dt, de/dt, di/dt, dOMEGA/dt, dw/dt et dM/dt respectivement
   
    eq1= (2/(n*a)) * dRdM    
    eq2= - ((1-e**2)*0.5)*dRdw/(n*e*a**2 ) + (1-e**2)*dRdM / (n*e*a**2)    
    eq3= - dRdOMEGA / (n*np.sin(i)*(a**2)*(1-e**2)**0.5) + np.cos(i)*dRdw/(n*np.sin(i)*(a**2)*(1-e**2)**0.5)   
    eq4= dRdi/(n*(a**2)*np.sin(i)*(1-e**2)**0.5)  
    eq5= ((1-e**2)**0.5)*dRde/(n*e*a**2) - np.cos(i)*dRdi/(n*(a**2)*np.sin(i)*(1-e**2)**0.5)
    eq6= n - 2*dRda/(n*a) - (1-e**2)*dRde/ (n*e*a**2)
    
    
    return np.array([eq1, eq2, eq3, eq4, eq5, eq6])


def rungekutta(  R,  t,  h):

    k1= h*dR( R , t )
    k2= h*dR( R+0.5*k1 , t+0.5*h )
    k3= h*dR( R+0.5*k2 , t+0.5*h )
    k4= h*dR( R+k3 , t+h )

    R= R + (1/6)*( k1+ 2*k2 + 2*k3 + k4 )
    return R
  """