# -*- coding: utf-8 -*-
"""
author -- Kazennova Daria
"""
from random import uniform
from math import exp, sqrt, pi, log, cos, sin
from matplotlib import pyplot
from pylab import polyfit

Diff_coeff = 0.2

# r - displacement, P - radius, t - time
def real_function(r, P, t):
    return exp(-(P*r**2)/(2*Diff_coeff*t))*sqrt(P/(2*Diff_coeff*pi*t))

def max_f(P, t):
    return sqrt(P/(2*Diff_coeff*pi*t))

def max_r(P, t):
    return sqrt((2*Diff_coeff*t*3*log(10))/P)

def random_pair(maxf, maxr):
    return (uniform(maxf/1000, maxf), uniform(0, maxr))
    
def correct_point(P, t):
    maxf = max_f(P, t)
    maxr = max_r(P, t)
    while True:
        point = random_pair(maxf, maxr)
        if real_function(point[1], P, t) >= point[0] :
            return point[1]
            break
        
def coordinates(r):
    phi = uniform(0, 2*pi)
    return (r*cos(phi), r*sin(phi))
    
def bounce(P, centerx, centery, centernewx, centernewy, R, flag):
    k = (centernewy - centery)/(centernewx - centerx)
    b = centery - centerx*k
    xplus1 = (-k*b + sqrt(R**2*k**2 + R**2 - b**2))/(1 + k**2)
    xminus1 = (-k*b - sqrt(R**2*k**2 + R**2 - b**2))/(1 + k**2)
    if centernewx < centerx :
        tempx = xminus1
    else :
        tempx = xplus1
    tempy = k*tempx + b
    cosinus = ((tempx - centerx)**2 + (tempy - centery)**2 + R**2 - centerx**2 -centery**2)/(2*sqrt((tempx - centerx)**2 + (tempy - centery)**2)*R)
    ct = R*cosinus - sqrt(R**2*cosinus**2 - 2*R*P + P**2)
    discr = (tempx - k*b + k*tempy)**2 - (1 + k**2)*(tempx**2 + b**2 - 2*b*tempy + tempy**2 - ct**2)
    xplus2 = (tempx - k*b + k*tempy + sqrt(discr))/(1 + k**2)
    xminus2 = (tempx - k*b + k*tempy - sqrt(discr))/(1 + k**2)
    if centerx < xplus2 < tempx or centerx > xplus2 > tempx :
        centertemp1x = xplus2
    else :
        centertemp1x = xminus2
    centertemp1y = k*centertemp1x + b
    k1 = centertemp1y/centertemp1x
    centertemp2x = (centertemp1y + centertemp1x/k1 - centernewy + k1*centernewx)/(k1 + 1/k1)
    centertemp2y = k1*centertemp2x + centernewy - k1*centernewx
    centerfinalx = 2*centertemp2x - centernewx
    centerfinaly = 2*centertemp2y - centernewy
    if sqrt(centerfinalx**2 + centerfinaly**2) > R - P :
        return bounce(P, centertemp1x, centertemp1y, centerfinalx, centerfinaly, R, flag)
    else :
        return [centerfinalx, centerfinaly, flag]    
        
def movement(P, centerx, centery, R):
    r = correct_point(P, 1)
    displacement = coordinates(r)
    centernewx = centerx + displacement[0]
    centernewy = centery + displacement[1]
    if sqrt(centernewx**2 + centernewy**2) > R - P :
        return bounce(P, centerx, centery, centernewx, centernewy, R, True)
    else :
        return [centernewx, centernewy, False]

def time_path_test(P, time, R):
    centernewx = [0]
    centernewy = [0]
    cn = [0,0,0]
    crossing = []
    for i in range(time+1):
        cn = movement(P, centernewx[i], centernewy[i], R)
        centernewx.append(cn[0])
        centernewy.append(cn[1])
        crossing.append(cn[2])
    r = [0]
    s = [0]
    for i in range(1, time+1):
        s.append(0)
        for j in range(1, time+2-i):
            s[i] += (centernewx[j-1]-centernewx[j-1+i])**2+(centernewy[j-1]-centernewy[j-1+i])**2
        r.append(s[i]/(time-i+1))
    all_crossing = sum(list(map(int, crossing)))
    log_t = [0]
    log_r = [0]
    for i in range(1, len(r)):
        log_t.append(log(i))
        log_r.append(log(r[i]))
    timemax = 15 + 1
    log_coeff = polyfit(log_t[1:timemax], log_r[1:timemax], 1)
    coeff = [log_coeff[0], exp(log_coeff[1])]
    coeff.append(all_crossing)
    return(coeff)
    
def average_coefficients(P, time, R):
    coeffs = [0, 0, 0]
    for i in range(1000):
        new_coeffs = time_path_test(P, time, R)
        coeffs[0] += new_coeffs[0]
        coeffs[1] += new_coeffs[1]
        coeffs[2] += new_coeffs[2]
    av_coeffs = [0, 0, 0]
    av_coeffs[0] = coeffs[0]/1000
    av_coeffs[1] = coeffs[1]/1000
    av_coeffs[2] = coeffs[2]/1000
    return(av_coeffs)

def R_test(P, time):
    R = [1.3]
    for i in range(1, 150):
        R.append(1.3 + i*0.1)
    print(R)        
    alpha_coeff = []
    a_coeff = []
    cross_coeff = []
    for i in range(150) :
        new_coeffs = average_coefficients(P, time, R[i])
        print(new_coeffs, R[i])
        alpha_coeff.append(new_coeffs[0])
        a_coeff.append(new_coeffs[1])
        cross_coeff.append(new_coeffs[2])
    pyplot.subplot(3, 1, 1)
    pyplot.plot(R, alpha_coeff, 'g.-')
    pyplot.title('Average coefficents (P = 1, t = 100)')
    pyplot.ylabel('alpha')
    pyplot.subplot(3, 1, 2)
    pyplot.plot(R, a_coeff, 'k.-')
    pyplot.ylabel('a')
    pyplot.subplot(3, 1, 3)
    pyplot.plot(R, cross_coeff, 'go')
    pyplot.xlabel('R')
    pyplot.ylabel('crossings')
    pyplot.savefig('av-coeff')

R_test(1, 100)
       