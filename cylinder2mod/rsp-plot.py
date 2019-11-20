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

    # Determine cylinder radius:
    A = REV * (2*EPS / (3 - EPS * EPSe)) ** (1. / 3.)
    # Determine half-length of the cylinder
    Htot = A/EPS
    He = A * EPSe
    Hc = Htot -He
    Hc = Hc


    # Parameters for nanorod cap
    aa = CAP**2
    bb = A**2

    for I in range(NGAUSS):
        CO = -X[I]
        SI = np.sqrt(1.-CO*CO)
        if Hc*SI > A*CO:
        # if np.abs(SI/CO) > np.abs(A/H):
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
            beta = -bb*2*Hc*CO
            gamma = bb*Hc**2 - aa*bb

            detsr = np.sqrt(beta**2 - (gamma)*(4*alpha))
            # RAD = (-beta - detsr)/(2*alpha)
            RAD = (-beta+np.sqrt(beta**2-4*alpha*gamma))/(2*alpha)
            # nom = -beta - detsr
            # psi = aa*SI*CO - bb*SI*CO
            # .replace('-beta - np.sqrt(beta**2 - (gamma)*(4*alpha))','nom') \
            # .replace('np.sqrt(beta**2 - (gamma)*(4*alpha))','detsr') \
            # valDR = -((2*alpha)*((-2*Hc*bb*SI - (-4*Hc**2*bb**2*SI*CO + (-gamma)*(8*psi)/2)/detsr)/(2*alpha) + (nom)*(-4*psi)/(2*alpha)**2)/(nom))
            # valDR = -(2*alpha)*((-2*Hc*bb*SI - (-4*Hc**2*bb**2*SI*CO + (-gamma)*(8*psi)/2)/detsr)/(2*alpha) + (nom)*(-4*psi)/(2*alpha)**2)/(nom)
            # valDR = -((2*aa*s2 + 2*bb*c2)*((-2*Hc*bb*SI - (-4*Hc**2*bb**2*SI*CO + (-Hc**2*bb + aa*bb)*(8*aa*SI*CO - 8*bb*SI*CO)/2)/np.sqrt(4*Hc**2*bb**2*c2 - (Hc**2*bb - aa*bb)*(4*aa*s2 + 4*bb*c2)))/(2*aa*s2 + 2*bb*c2) + (2*Hc*bb*CO - np.sqrt(4*Hc**2*bb**2*c2 - (Hc**2*bb - aa*bb)*(4*aa*s2 + 4*bb*c2)))*(-4*aa*SI*CO + 4*bb*SI*CO)/(2*aa*s2 + 2*bb*c2)**2)/(2*Hc*bb*CO - np.sqrt(4*Hc**2*bb**2*c2 - (Hc**2*bb - aa*bb)*(4*aa*s2 + 4*bb*c2)))
            #           )
            valDR = (2*aa*s2 + 2*bb*c2)*(-(-2*Hc*bb*SI - (-4*Hc**2*bb**2*SI*CO + (-Hc**2*bb + aa*bb)*(8*aa*SI*CO - 8*bb*SI*CO)/2)/np.sqrt(4*Hc**2*bb**2*c2 - (Hc**2*bb - aa*bb)*(4*aa*s2 + 4*bb*c2)))/(2*aa*s2 + 2*bb*c2) - (2*Hc*bb*CO - np.sqrt(4*Hc**2*bb**2*c2 - (Hc**2*bb - aa*bb)*(4*aa*s2 + 4*bb*c2)))*(-4*aa*SI*CO + 4*bb*SI*CO)/(2*aa*s2 + 2*bb*c2)**2)/(2*Hc*bb*CO - np.sqrt(4*Hc**2*bb**2*c2 - (Hc**2*bb - aa*bb)*(4*aa*s2 + 4*bb*c2)))

            # valDR = -(
            #            (-2*Hc*bb*SI - (-4*Hc**2*bb**2*SI*CO -4*gamma*psi
            #                         )
            #                         /detsr
            #            )
            #            - 2*nom*psi/alpha
            #          )/nom

        R[I] = RAD**2
        R[NG-I-1] = R[I]          # using mirror symmetry
        DR[I] = valDR
        DR[NG-I-1] = -DR[I]       # using mirror symmetry

    return R, DR

def rsp_nanorod2(X, NG, NGAUSS, REV, EPS, EPSe):
    R = np.zeros(NG)
    DR = np.zeros(NG)
    # Determine cylinder radius:
    A = REV * (2*EPS / (3 - EPS * EPSe)) ** (1. / 3.)
    # Determine half-length of the cylinder
    H = A/EPS
    He = A * EPSe
    Hc = H -He
    EPSc = Hc/A

    for I in range(NGAUSS):
        CO = -X[I]
        SI = np.sqrt(1.-CO*CO)
        # if np.abs(CO/SI) < np.abs(Hc/A):
        if Hc*SI > A*CO:
            # Along the circular surface:
            RAD = A/SI
            RTHET = -A*CO/(SI*SI)
            valDR = -RTHET/RAD

        else:
            # Along elliptic cap
            c2 = CO**2
            s2 = SI**2
            # Solution of square euation of ellipse move from the origin
            alpha = np.sqrt(
                            (EPSe**2-EPSc**2)*s2
                            +
                            c2
                           )
            beta = EPSe**2*s2 + c2
            RAD = (Hc*CO +  #
                     He*alpha)/beta
            valDR = (-alpha*Hc*SI
                        + He*(EPSe**2-EPSc**2-1)*SI*CO
                    )/(
                        He*alpha**2 +  #
                        alpha*Hc*CO
                    )-(
                        2*(EPSe**2-1.)*SI*CO/beta
                      )
            valDR *= -1



        R[I] = RAD**2
        R[NG-I-1] = R[I]          # using mirror symmetry
        DR[I] = valDR
        DR[NG-I-1] = -DR[I]       # using mirror symmetry

    return R, DR



NGAUSS = 200
NG = 2*NGAUSS
REV = 2.  # Equivalent radius
EPS = 0.5
# EPSe = 1/EPS
EPSe = 1.1


A1 = REV * (2*EPS / (3 - EPS * EPSe)) ** (1. / 3.)
H1 = A1/EPS
CAP = A1 * EPSe

theta = np.linspace(0, np.pi, NG)[::-1]
theta2 = np.linspace(0, np.pi/10, NG)[::-1]
# theta = np.linspace(-np.pi/2., np.pi/2., NG)
X = np.cos(theta)
R1, DR1 = rsp1(X, NG, NGAUSS, REV, EPS)
R3, DR3 = rsp3(X, NG, NGAUSS, REV, EPS)
R, DR = rsp_nanorod(X, NG, NGAUSS, REV, EPS, CAP)
Rc, DRc = rsp_nanorod2(X, NG, NGAUSS, REV, EPS, EPSe)

# isDecart = True
isDecart = False
if isDecart:
    xxa = R**(1./2.)*np.cos(theta)
    yya = R**(1./2.)*np.sin(theta)
    xxac = Rc**(1./2.)*np.cos(theta)
    yyac = Rc**(1./2.)*np.sin(theta)
    xxa1 = R1**(1./2.)*np.cos(theta)
    yya1 = R1**(1./2.)*np.sin(theta)
    xxa3 = R3**(1./2.)*np.cos(theta)
    yya3 = R3**(1./2.)*np.sin(theta)
    xcirc = H1 * np.sin(theta2)
    ycirc = H1 * np.cos(theta2)

    xx3 = DR3*np.cos(theta)*R3**(1./2.)
    yy3 = DR3*np.sin(theta)*R3**(1./2.)
    xx1 = DR1*np.cos(theta)*R1**(1./2.)
    yy1 = DR1*np.sin(theta)*R1**(1./2.)
    xx = DR*np.cos(theta)*R**(1./2.)
    yy = DR*np.sin(theta)*R**(1./2.)
    xxc = DRc*np.cos(theta)*Rc**(1./2.)
    yyc = DRc*np.sin(theta)*Rc**(1./2.)
else:
    xxa = theta
    xxac = theta
    xxa1 = theta
    xxa3 = theta
    xcirc = theta

    yya = R**(1./2.)
    yyac = Rc**(1./2.)
    yya1 = R1**(1./2.)
    yya3 = R3**(1./2.)
    ycirc = H1*np.ones(len(theta))



    xx = theta
    xxc = theta
    xx1 = theta
    xx3 = theta

    yy = DR*R**(1./2.)
    yyc = DRc*Rc**(1./2.)
    yy1 = DR1*R1**(1./2.)
    yy3 = DR3*R3**(1./2.)


color = np.linspace(0, 1, len(xx))

# if True:
# if False:
fig, ax = plt.subplots(1,2)
ax[0].scatter(xxc, yyc, c=color, s=80, cmap='Set1')
ax[0].scatter(xx, yy, c=color, s=20, cmap='tab10', marker='x')
ax[0].scatter(xx3, yy3, c=color, s=10, cmap='jet')
ax[0].scatter(xx1, yy1, c=color, s=1, cmap='plasma_r')
    # ax[0].set_aspect(1.0)
    # ax[1].scatter(xx3, yy3, c=color, cmap='gray')
    # ax[1].set_aspect(1.0)
    # ax[2].scatter(xx1, yy1, c=color, cmap='gray')
    # ax[2].set_aspect(1.0)
    # for i in range(3):
    #     ax[i].set_xlim(-0.33,1.1)
    #     ax[i].set_ylim(-2,2)
# else:
#     fig, ax = plt.subplots(1,1)
ax[1].plot(xcirc,ycirc)
ax[1].scatter(xxac, yyac, c=color, s=80,cmap='Set1')
ax[1].scatter(xxa, yya, c=color, s=20,cmap='tab10', marker ='x')
ax[1].scatter(xxa1, yya1, c=color,s=1, cmap='plasma_r')
ax[1].scatter(xxa3, yya3, c=color,s=10, cmap='jet_r')
# ax[0].set_aspect(1.0)
# ax[1].scatter(xxa3, yya3, c=color, cmap='gray')
ax[1].set_aspect(1.0)
# ax[2].scatter(xxa1, yya1, c=color, cmap='gray')
# ax[2].set_aspect(1.0)

import sympy as sp

x, aa, bb, Hc = sp.symbols('x aa bb Hc', real=True)
c2 = sp.cos(x)**2
s2 = sp.sin(x)**2
CO = sp.cos(x)
SI = sp.sin(x)

alpha = bb*c2 + aa*s2
beta = -2*Hc*bb*CO
gamma = Hc**2*bb - aa*bb
RAD = (-beta-sp.sqrt(beta**2-4*alpha*gamma))/(2*alpha)

rdr = str(sp.diff(-RAD,x)/RAD) \
    .replace('sin(x)**2','s2') \
    .replace('cos(x)**2','c2') \
    .replace('sin(x)','SI') \
    .replace('cos(x)','CO') \
    .replace('sqrt', 'np.sqrt') \
    # .replace('2*aa*s2 + 2*bb*c2','2*alpha') \
    # .replace('4*aa*s2 + 4*bb*c2','4*alpha') \
    # .replace('2*Hc*bb*CO','-beta')\
    # .replace('4*Hc**2*bb**2*c2','beta**2')\
    # .replace('-Hc**2*bb + aa*bb','-gamma') \
    # .replace('Hc**2*bb - aa*bb','gamma')\
    # .replace('-beta - np.sqrt(beta**2 - (gamma)*(4*alpha))','nom')\
    # .replace('np.sqrt(beta**2 - (gamma)*(4*alpha))','detsr')\
    # .replace('8*aa*SI*CO - 8*bb*SI*CO', '8*psi')\
    # .replace('-4*aa*SI*CO + 4*bb*SI*CO', '-4*psi')


print(rdr)


plt.show()
