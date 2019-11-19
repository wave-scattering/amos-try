#!/usr/bin/python3
# -*- coding: UTF-8 -*-
from matplotlib import pyplot as plt
import numpy as np
import math as m

def rsp1 (X,NG,NGAUSS,REV,EPS):
    R = np.zeros(NG)
    DR = np.zeros(NG)
    A = REV*EPS**(1./3.)
    AA = A*A
    EE = EPS*EPS
    EE1 = EE-1.

    for I in range(NGAUSS):
        C = X[I]
        CC = C*C
        SS = 1.-CC
        S = np.sqrt(SS)                 #  = \sin\theta
        RR = 1./(SS+EE*CC)
        R[I] = AA*RR
        R[NG-I-1] = R[I]
        DR[I] = RR*C*S*(EE-1)
        DR[NG-I-1] = -DR[I]
    return R, DR


def rsp3(X, NG, NGAUSS, REV, EPS):
    R = np.zeros(NG)
    DR = np.zeros(NG)
    # Determine half-length of the cylinder
    H = REV*( (2./(3.*EPS*EPS))**(1./3.) )
    # Determine cylinder radius:
    A = H*EPS
    for I in range(NGAUSS):
        CO = -X[I]
        SI = np.sqrt(1.-CO*CO)
        if np.abs(SI/CO) > np.abs(A/H):
            # Along the circular surface:
            RAD = A/SI
            RTHET = -A*CO/(SI*SI)
        else:
            # Along the plane cuts:
            RAD = H/CO
            RTHET = H*SI/(CO*CO)
        R[I] = RAD*RAD
        R[NG-I-1] = R[I]          # using mirror symmetry
        DR[I] = -RTHET/RAD
        DR[NG-I-1] = -DR[I]       # using mirror symmetry
    return R, DR


def rsp_nanorod(X, NG, NGAUSS, REV, EPS, CAP):
    R = np.zeros(NG)
    DR = np.zeros(NG)
    # Determine half-length of the cylinder
    H = REV*( (2./(3.*EPS*EPS))**(1./3.) )
    # Determine cylinder radius:
    A = H*EPS

    if isOnlySpheroid:
        A = REV*EPS**(1./3.)
        H = 0.0001

    # Parameters for nanorod cap
    aa = CAP**2
    bb = A**2

    for I in range(NGAUSS):
        CO = -X[I]
        SI = np.sqrt(1.-CO*CO)
        if np.abs(SI/CO) > np.abs(A/H):
            # Along the circular surface:
            RAD = A/SI
            RTHET = -A*CO/(SI*SI)
            valDR = -RTHET/RAD

        else:
            # Along elliptic cap
            c2 = CO**2
            s2 = SI**2
            # Solution of square euation of ellipse move from the origin
            alpha = bb*c2 + aa*s2
            beta = -bb*2*H*CO
            gamma = bb*H**2 - aa*bb
            detsr = np.sqrt(beta**2 - (gamma)*(4*alpha))
            nom = -beta - detsr
            RAD = nom/(2*alpha)
            psi = aa*SI*CO - bb*SI*CO

            valDR = -(
                       (-2*H*bb*SI - (-4*H**2*bb**2*SI*CO -4*gamma*psi
                                    )
                                    /detsr
                       )
                       - 2*nom*psi/alpha
                     )/nom

        R[I] = RAD**2
        R[NG-I-1] = R[I]          # using mirror symmetry
        DR[I] = valDR
        DR[NG-I-1] = -DR[I]       # using mirror symmetry

    return R, DR


NGAUSS = 100
NG = 2*NGAUSS
REV = 2.  # Equivalent radius
EPS = 1
CAP = 0.0004
isOnlySpheroid = True
isOnlySpheroid = False
if isOnlySpheroid:
    CAP = REV*EPS**(1./3.)/2

theta = np.linspace(0, np.pi, NG)[::-1]
# theta = np.linspace(-np.pi/2., np.pi/2., NG)
X = np.cos(theta)
R1, DR1 = rsp1(X, NG, NGAUSS, REV, EPS)
R3, DR3 = rsp3(X, NG, NGAUSS, REV, EPS)
R, DR = rsp_nanorod(X, NG, NGAUSS, REV, EPS, CAP)
print(DR)

xxa = R**(1./2.)*np.cos(theta)
yya = R**(1./2.)*np.sin(theta)
xxa1 = R1**(1./2.)*np.cos(theta)
yya1 = R1**(1./2.)*np.sin(theta)
xxa3 = R3**(1./2.)*np.cos(theta)
yya3 = R3**(1./2.)*np.sin(theta)


xx3 = DR3*np.cos(theta)*R3**(1./2.)
yy3 = DR3*np.sin(theta)*R3**(1./2.)
xx1 = DR1*np.cos(theta)*R1**(1./2.)
yy1 = DR1*np.sin(theta)*R1**(1./2.)
xx = DR*np.cos(theta)*R**(1./2.)
yy = DR*np.sin(theta)*R**(1./2.)

color = np.linspace(0, 1, len(xx))

if True:
# if False:
    fig, ax = plt.subplots(1,1)
    ax.scatter(xx3, yy3, c=color, s=100, cmap='jet')
    ax.scatter(xx1, yy1, c=color, s=100, cmap='plasma_r')
    ax.scatter(xx, yy, c=color, s=20, cmap='gray')
    # ax[0].set_aspect(1.0)
    # ax[1].scatter(xx3, yy3, c=color, cmap='gray')
    # ax[1].set_aspect(1.0)
    # ax[2].scatter(xx1, yy1, c=color, cmap='gray')
    # ax[2].set_aspect(1.0)
    # for i in range(3):
    #     ax[i].set_xlim(-0.33,1.1)
    #     ax[i].set_ylim(-2,2)
else:
    fig, ax = plt.subplots(1,1)
    ax.scatter(xxa, yya, c=color, s=80,cmap='gray')
    ax.scatter(xxa1, yya1, c=color,s=20, cmap='plasma_r')
    ax.scatter(xxa3, yya3, c=color,s=20, cmap='jet_r')
    # ax[0].set_aspect(1.0)
    # ax[1].scatter(xxa3, yya3, c=color, cmap='gray')
    # ax[1].set_aspect(1.0)
    # ax[2].scatter(xxa1, yya1, c=color, cmap='gray')
    # ax[2].set_aspect(1.0)

import sympy as sp

x, aa, bb, H = sp.symbols('x aa bb H', real=True)
c2 = sp.cos(x)**2
s2 = sp.sin(x)**2
CO = sp.cos(x)
SI = sp.sin(x)

alpha = bb*c2 + aa*s2
beta = -2*H*bb*CO
gamma = H**2*bb - aa*bb
RAD = (-beta-sp.sqrt(beta**2-4*alpha*gamma))/(2*alpha)

rdr = str(sp.diff(RAD,x)/RAD) \
    .replace('sin(x)**2','s2') \
    .replace('cos(x)**2','c2') \
    .replace('sin(x)','SI') \
    .replace('cos(x)','CO') \
    .replace('sqrt', 'np.sqrt') \
    .replace('2*aa*s2 + 2*bb*c2','2*alpha') \
    .replace('4*aa*s2 + 4*bb*c2','4*alpha') \
    .replace('2*H*bb*CO','-beta')\
    .replace('4*H**2*bb**2*c2','beta**2')\
    .replace('-H**2*bb + aa*bb','-gamma') \
    .replace('H**2*bb - aa*bb','gamma')\
    .replace('-beta - np.sqrt(beta**2 - (gamma)*(4*alpha))','nom')\
    .replace('np.sqrt(beta**2 - (gamma)*(4*alpha))','detsr')\
    .replace('8*aa*SI*CO - 8*bb*SI*CO', '8*psi')\
    .replace('-4*aa*SI*CO + 4*bb*SI*CO', '-4*psi')


#print(rdr)


plt.show()
