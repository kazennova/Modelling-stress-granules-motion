# -*- coding: utf-8 -*-
"""
author -- Kazennova Daria
"""
from random import uniform
from math import exp, sqrt, pi, log, cos, sin
from pylab import  *

Diff_coeff = 1.5
Domens = []
big_domens = 5
small_domens = 30
def create_domen(r_min, r_max):
    r = uniform(r_min, r_max)
    phi = uniform(0, 2*pi)
    center = uniform(0, 20)
    return [center*cos(phi), center*sin(phi), r]

def create_big_domens():
    new_domen = []
    new_domen = create_domen(3, 4)
    Domens.append(new_domen)
    for i in range(big_domens - 1):
        while True:
            k = 0
            new_domen = create_domen(3, 4)
            for j in range(len(Domens)):
                if sqrt((new_domen[0]-Domens[j][0])**2 + (new_domen[1]-Domens[j][1])**2) > max(new_domen[2], Domens[j][2]):
                    k = k + 1
            if k == len(Domens) :
                Domens.append(new_domen)
                break

def create_small_domens():
    new_domen = []
    for i in range(small_domens):
        while True:
            k = 0
            new_domen = create_domen(0.2, 0.5)
            for j in range(len(Domens)):
                if sqrt((new_domen[0]-Domens[j][0])**2 + (new_domen[1]-Domens[j][1])**2) > max(new_domen[2], Domens[j][2]):
                    k = k + 1
            if k == len(Domens) :
                Domens.append(new_domen)
                break

def create_domens():
    Domens.clear()
    create_big_domens()
    create_small_domens()
    return(Domens)

def not_available_point(P, centerx, centery, centernewx, centernewy, R):
    if sqrt(centernewx**2 + centernewy**2) <= R - P :
        for i in range(len(Domens)) :
            if distance_to_domen(P, centerx, centery, centernewx, centernewy, Domens[i][2], Domens[i][0], Domens[i][1], R) < 2*R :
                return(True)
        return(False)
    else :
        return(True)

def start_point(P, R):
    r = uniform(0, R - P)
    phi = uniform(0, 2*pi)
    point = [r*cos(phi), r*sin(phi)]
    if sqrt(point[0]**2 + point[1]**2) <= R - P :
        for i in range(len(Domens)) :
            if P + Domens[i][2] > sqrt((point[0] - Domens[i][0])**2 + (point[1] - Domens[i][1])**2) :
                return(start_point(P, R))
        return(point)
    else :
        return(start_point(P, R))

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
    
def distance_to_domen(P, centerx, centery, centernewx, centernewy, r, domenx, domeny, R):
    k = (centernewy - centery)/(centernewx - centerx)
    b = -k*centerx + centery
    k1 = -1/k
    b1 = domeny - domenx*k1
    prx = (b1-b)/(k-k1)
    pry = (k*b1-k1*b)/(k-k1)
    dtd = 2*R
    if sqrt((domenx - prx)**2 + (domeny - pry)**2) < (P + r) :
        dist = sqrt((P+r)**2-(domenx-prx)**2-(domeny-pry)**2)
        if centernewx < centerx :
            tempx = prx + dist/sqrt(1+k**2)  
        else :
            tempx = prx - dist/sqrt(1+k**2)
        tempy = k*tempx + b
        if centerx < tempx < centernewx or centerx > tempx > centernewx :
            dtd = sqrt((tempx - centerx)**2 + (tempy - centery)**2)
    return(dtd)
    
def bounce_from_domen(P, centerx, centery, centernewx, centernewy, r, domenx, domeny):
    k = (centernewy - centery)/(centernewx - centerx)
    b = -k*centerx + centery
    k1 = -1/k
    b1 = domeny - domenx*k1
    prx = (b1-b)/(k-k1)
    pry = (k*b1-k1*b)/(k-k1)
    dist = sqrt((P+r)**2-(domenx-prx)**2-(domeny-pry)**2)
    if centernewx < centerx :
        tempx = prx + dist/sqrt(1+k**2)  
    else :
        tempx = prx - dist/sqrt(1+k**2)   
    tempy = k*tempx + b
    k2 = (tempy - domeny)/(tempx - domenx)
    b2 = centernewy - k2*centernewx
    k3 = -1/k2
    b3 = tempy - k3*tempx
    centerfinalx = 2*(b3 - b2)/(k2 - k3) - centernewx
    centerfinaly = 2*(k2*b3 - k3*b2)/(k2 - k3) - centernewy
    return [tempx, tempy, centerfinalx, centerfinaly]
   
def distance_to_bound(P, centerx, centery, centernewx, centernewy, R) :
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
    dtb = sqrt((centerx - centertemp1x)**2 + (centery - centertemp1y)**2)
    return(dtb)

def bounce_from_bound(P, centerx, centery, centernewx, centernewy, R) :
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
    return [centertemp1x, centertemp1y, centerfinalx, centerfinaly]    
   
def movement(P, centerx, centery, R):
    r = correct_point(P, 1.59)
    displacement = coordinates(r)
    centernewx = centerx + displacement[0]
    centernewy = centery + displacement[1]
    while not_available_point(P, centerx, centery, centernewx, centernewy, R) :
        distances = []
        for k in range(len(Domens)) :
            distances.append(distance_to_domen(P, centerx, centery, centernewx, centernewy, Domens[k][2], Domens[k][0], Domens[k][1], R))
        distances.append(distance_to_bound(P, centerx, centery, centernewx, centernewy, R))
        index_min_distance = distances.index(sorted(distances)[0])
        if index_min_distance == len(distances) - 1 :
            centers = bounce_from_bound(P, centerx, centery, centernewx, centernewy, R)
        else :
            centers = bounce_from_domen(P, centerx, centery, centernewx, centernewy, Domens[index_min_distance][2], Domens[index_min_distance][0], Domens[index_min_distance][1])
        centerx = centers[0]
        centery = centers[1]
        centernewx = centers[2]
        centernewy = centers[3]
    return [centernewx, centernewy]

def time_path_test(P, time, R):
    centernew = start_point(P, R)
    centernewx = [centernew[0]]
    centernewy = [centernew[1]]
    cn = [0,0]
    for i in range(time+1):
        cn = movement(P, centernewx[i], centernewy[i], R)
        centernewx.append(cn[0])
        centernewy.append(cn[1])
    r = [0]
    s = [0]
    for i in range(1, time+1):
        s.append(0)
        for j in range(1, time+2-i):
            s[i] += (centernewx[j-1]-centernewx[j-1+i])**2+(centernewy[j-1]-centernewy[j-1+i])**2
        r.append(s[i]/(time-i+1))
    log_t = [0]
    log_r = [0]
    for i in range(1, len(r)):
        log_t.append(log(i))
        log_r.append(log(r[i]))
    timemax15 = 15 + 1
    log_coeff15 = polyfit(log_t[1:timemax15], log_r[1:timemax15], 1)
    coeff15 = [log_coeff15[0], exp(log_coeff15[1])]
    return(coeff15)

def P_test(time, R, dup):
    P = [0.149,0.252,0.259,0.299,0.329,0.329,0.329,0.334,0.334,0.343,0.357,0.366,0.370,0.370,0.370,0.374,0.399,0.407,0.411,0.411,0.418,0.418,0.433,0.433,0.444,0.455,0.465,0.465,0.465,0.472,0.479,0.492,0.498,0.505,0.508,0.511,0.514,0.514,0.520,0.532,0.538,0.573,0.578,0.597,0.605,0.610,0.621,0.623,0.631,0.643,0.653,0.660,0.663,0.693,0.707,0.714,0.723,0.753,0.765,0.790,0.806,0.816,0.816,0.820,0.827,0.827,0.835,0.839,0.843,0.856,0.863,0.872,0.881,0.915,0.918,0.927,0.934,0.939,0.944,0.946,0.956,0.977,0.980,0.980,0.997,0.997,1.037,1.042,1.056,1.069,1.076,1.079,1.095,1.116,1.123,1.138,1.142,1.163,1.205,1.221,1.235,1.239,1.251,1.270,1.288,1.289,1.305,1.332,1.333,1.374,1.395,1.437,1.444,1.454,1.456,1.542,1.757,2.173]
    P = P*dup
    P.sort()    
    alpha_coeff = [0]*len(P)
    a_coeff = [0]*len(P)
    for i in range(len(P)) :
        create_domens()
        new_coeffs = time_path_test(P[i], time, R)
        alpha_coeff[i] = new_coeffs[0]
        a_coeff[i] = new_coeffs[1]
    bin1 = 25*dup
    bin2 = bin1 + 25*dup
    bin3 = bin2 + 23*dup
    bin4 = bin3 + 24*dup
    number = [0]*15
    for i in range(len(P)) :
        if i < bin1 :
            if a_coeff[i] <= 0.02 :
                number[0] += 1
            else :
                if a_coeff[i] <= 0.06 :
                    number[1] +=1
                else :
                    number[2] +=1
        if bin1 <= i < bin2 :
            if a_coeff[i] <= 0.02 :
                number[3] += 1
            else :
                if a_coeff[i] <= 0.06 :
                    number[4] +=1
                else :
                    number[5] +=1      
        if bin2 <= i < bin3 :
            if a_coeff[i] <= 0.02 :
                number[6] += 1
            else :
                if a_coeff[i] <= 0.06 :
                    number[7] +=1
                else :
                    number[8] +=1   
        if bin3 <= i < bin4 :
            if a_coeff[i] <= 0.02 :
                number[9] += 1
            else :
                if a_coeff[i] <= 0.06 :
                    number[10] +=1
                else :
                    number[11] +=1        
        if i >= bin4 :
            if a_coeff[i] <= 0.02 :
                number[12] += 1
            else :
                if a_coeff[i] <= 0.06 :
                    number[13] +=1
                else :
                    number[14] +=1 
    hist = [0]*15
    for k in range(len(number)) :
        if k < 3 :
            hist[k] = number[k]/sum(number[:3])
        if 3 <= k < 6 :
            hist[k] = number[k]/sum(number[3:6])
        if 6 <= k < 9 :
            hist[k] = number[k]/sum(number[6:9]) 
        if 9 <= k < 12 :
            hist[k] = number[k]/sum(number[9:12])
        if  k >= 12 :
            hist[k] = number[k]/sum(number[12:]) 
    hist_real = [0.4583333333333333, 0.3333333333333333, 0.20833333333333334, 0.75, 0.125, 0.125, 0.625, 0.3333333333333333, 0.041666666666666664, 0.6956521739130435, 0.08695652173913043, 0.21739130434782608, 0.8260869565217391, 0.13043478260869565, 0.043478260869565216]
    sigma = 0
    for j in range(15) :
        sigma += (hist[j]-hist_real[j])**2
    sigma = sqrt(sigma/15)
    print(3, 4, big_domens, Diff_coeff, sigma, hist)

for big_domens in [5, 7] :
    for Diff_coeff in [0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08] :
        P_test(100, 20, 20)
