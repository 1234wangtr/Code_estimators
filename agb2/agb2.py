
####################      Readme #####################################
'''
You can run our tool on Windows using the command line interface. Our script runs with Python 3 and the numpy library. The input format for our tool includes two types: one for regular LPN over F2 and another for LPN over a larger field, where |F| > 2. Assuming ``AGBscript.py'' (our tool) is in the "C:" directory.

==================Parameters==========
The order of command line parameters (n, k, t, q) is fixed, and we will explain their meanings.
n: the length of noise
k: dimension of LPN
t: Hamming weight of a noise vector
q: for regular LPN over a larger field

================ input format for for regular LPN over F2:
python C:\\AGBscript.py n=1024 k=652 t=57

================ input format for for regular LPN over F2:
python C:\\AGBscript.py n=1024 k=652 t=57 q

'''

import math
import numpy as np
import decimal
import time
import datetime

import sys
import re
from decimal import Decimal




repeat = 10
search = 1000
decimal.getcontext().prec = 170
log10Two = decimal.Decimal(2).log10()


#####################      The recent algebraic attack (AGB),  in EC23      ###########################

####################### The cost of AGB over F_q

# compute degree_conjforq = d_{wit,(f,mu)} for larger field, d = 2,3,4,1000000

def degree_conjforq(n, k, h, f, mu):
    n = decimal.Decimal(n)
    k = decimal.Decimal(k)
    h = decimal.Decimal(h)
    f = decimal.Decimal(f)
    mu = decimal.Decimal(mu)
    beta = decimal.Decimal(math.floor(n / h))

    a1 = -(n - k - h)
    a2 = (n - k - h) * (n - k - h - 1) / 2
    a3 = -(n - k - h) * (n - k - h - 1) * (n - k - h - 2) / 6
    a4 = (n - k - h) * (n - k - h - 1) * (n - k - h - 2) * (n - k - h - 3) / 24

    b1 = (beta - mu - 1) * f
    b2 = (beta - mu - 1) ** 2 * f * (f - 1) / 2
    b3 = (beta - mu - 1) ** 3 * f * (f - 1) * (f - 2) / 6
    b4 = (beta - mu - 1) ** 4 * f * (f - 1) * (f - 2) * (f - 3) / 24

    c1 = (beta - 1) * (h - f)
    c2 = (beta - 1) ** 2 * (h - f) * (h - f - 1) / 2
    c3 = (beta - 1) ** 3 * (h - f) * (h - f - 1) * (h - f - 2) / 6
    c4 = (beta - 1) ** 4 * (h - f) * (h - f - 1) * (h - f - 2) * (h - f - 3) / 24

    d2 = a1 * b1 + a1 * c1 + b1 * c1 + a2 + b2 + c2
    d3 = c3 + b1 * c2 + b2 * c1 + b3 + a1 * (b1 * c1 + b2 + c2) + a2 * (b1 + c1) + a3
    d4 = c4 + b1 * c3 + b2 * c2 + b3 * c1 + b4 + a1 * (b1 * c2 + b2 * c1 + b3 + c3) + a2 * (b2 + c2 + b1 * c1) + a3 * (
            b1 + c1) + a4

    if d2 < 1:
        return 2

    if d3 + d2 < 1:
        return 3

    if d4 + d3 + d2 < 1:
        return 4

    return 30


# compute M for parameter f and mu

def subAGBforq(n, k, h, f, mu):
    n = decimal.Decimal(n)
    k = decimal.Decimal(k)
    h = decimal.Decimal(h)
    f = decimal.Decimal(f)
    mu = decimal.Decimal(mu)
    beta = decimal.Decimal(math.floor(n / h))

    a1 = (beta - mu - 1) * f
    a2 = (beta - mu - 1) ** 2 * f * (f - 1) / 2
    a3 = (beta - mu - 1) ** 3 * f * (f - 1) * (f - 2) / 6
    a4 = (beta - mu - 1) ** 4 * f * (f - 1) * (f - 2) * (f - 3) / 24

    b1 = (beta - 1) * (h - f)
    b2 = (beta - 1) ** 2 * (h - f) * (h - f - 1) / 2
    b3 = (beta - 1) ** 3 * (h - f) * (h - f - 1) * (h - f - 2) / 6
    b4 = (beta - 1) ** 4 * (h - f) * (h - f - 1) * (h - f - 2) * (h - f - 3) / 24

    c1 = (h - 1)
    c2 = (h) * (h - 1) / 2
    c3 = (h + 1) * (h) * (h - 1) / 6
    c4 = (h + 2) * (h + 1) * (h) * (h - 1) / 24

    d = degree_conjforq(n, k, h, f, mu)
    d2 = a1 + b1 + c1 + a1 * b1 + a1 * c1 + b1 * c1 + a2 + b2 + c2
    if d == 2:
        M = d2

    elif d == 3:
        d3 = c3 + b1 * c2 + b2 * c1 + b3 + a1 * (b1 * c1 + b2 + c2) + a2 * (b1 + c1) + a3
        M = d2 + d3

    elif d == 4:
        d3 = c3 + b1 * c2 + b2 * c1 + b3 + a1 * (b1 * c1 + b2 + c2) + a2 * (b1 + c1) + a3
        d4 = c4 + b1 * c3 + b2 * c2 + b3 * c1 + b4 + a1 * (b1 * c2 + b2 * c1 + b3 + c3) + a2 * (
                b2 + c2 + b1 * c1) + a3 * (b1 + c1) + a4
        M = d2 + d3 + d4

    else:
        return decimal.Decimal(2) ** 300

    return 2 * M.log10() / log10Two + (3 * (k + 1 - f * mu)).log10() / log10Two - f * (
            1 - mu / beta).log10() / log10Two


# The cost of AGB over F_q
def AGBforq(n, k, h):
    n = decimal.Decimal(n)
    k = decimal.Decimal(k)
    finalcost = 999999
    finalf = 99999
    finalmu = 999999
    for f in range(0, h + 1):
        for mu in range(0, math.floor(n / h) + 1):
            if f * mu < k + 1:
                subcost = subAGBforq(n, k, h, f, mu)
                if finalcost > subcost:
                    finalcost = subcost
                    finalf = f
                    finalmu = mu
    # print(degree_conjforq(n, k, h, finalf, finalmu))
    # print(finalf)
    # print(finalmu)
    # print("   for larger field   n = " + str(n) + "     k = " + str(k) + "    h = " + str(h) + "  cost  " + str(    round(finalcost, 8)))
    return round(finalcost, 2)


####################### The cost of AGB over F_2

# compute degree_conjforq = d_{wit,(f,mu)} for F2, d = 2,3,4,1000000
def degree_conjfor2(n, k, h, f, mu):
    n = decimal.Decimal(n)
    k = decimal.Decimal(k)
    h = decimal.Decimal(h)
    f = decimal.Decimal(f)
    mu = decimal.Decimal(mu)
    beta = decimal.Decimal(math.floor(n / h))
    a1 = (beta - mu - 1) * f
    a2 = (beta - mu - 1) ** 2 * f * (f - 1) / 2
    a3 = (beta - mu - 1) ** 3 * f * (f - 1) * (f - 2) / 6
    a4 = (beta - mu - 1) ** 4 * f * (f - 1) * (f - 2) * (f - 3) / 24

    b1 = (beta - 1) * (h - f)
    b2 = (beta - 1) ** 2 * (h - f) * (h - f - 1) / 2
    b3 = (beta - 1) ** 3 * (h - f) * (h - f - 1) * (h - f - 2) / 6
    b4 = (beta - 1) ** 4 * (h - f) * (h - f - 1) * (h - f - 2) * (h - f - 3) / 24

    c1 = -(n - k - 1)
    c2 = (n - k) * (n - k - 1) / 2
    c3 = -(n - k + 1) * (n - k) * (n - k - 1) / 6
    c4 = (n - k + 2) * (n - k + 1) * (n - k) * (n - k - 1) / 24

    d2 = a1 + b1 + c1 + a1 * b1 + a1 * c1 + b1 * c1 + a2 + b2 + c2
    d3 = c3 + b1 * c2 + b2 * c1 + b3 + a1 * (b1 * c1 + b2 + c2) + a2 * (b1 + c1) + a3
    d4 = c4 + b1 * c3 + b2 * c2 + b3 * c1 + b4 + a1 * (b1 * c2 + b2 * c1 + b3 + c3) + a2 * (b2 + c2 + b1 * c1) + a3 * (
            b1 + c1) + a4

    if d2 < 1:
        return 2

    if d3 + d2 < 1:
        return 3

    if d4 + d3 + d2 < 1:
        return 4

    return 30


# compute M for parameter f and mu
def subAGBfor2(n, k, h, f, mu):
    n = decimal.Decimal(n)
    k = decimal.Decimal(k)
    h = decimal.Decimal(h)
    f = decimal.Decimal(f)
    mu = decimal.Decimal(mu)
    beta = decimal.Decimal(math.floor(n / h))

    a1 = (beta - mu - 1) * f
    a2 = (beta - mu - 1) ** 2 * f * (f - 1) / 2
    a3 = (beta - mu - 1) ** 3 * f * (f - 1) * (f - 2) / 6
    a4 = (beta - mu - 1) ** 4 * f * (f - 1) * (f - 2) * (f - 3) / 24

    b1 = (beta - 1) * (h - f)
    b2 = (beta - 1) ** 2 * (h - f) * (h - f - 1) / 2
    b3 = (beta - 1) ** 3 * (h - f) * (h - f - 1) * (h - f - 2) / 6
    b4 = (beta - 1) ** 4 * (h - f) * (h - f - 1) * (h - f - 2) * (h - f - 3) / 24

    d = degree_conjfor2(n, k, h, f, mu)
    d2 = a1 + b1 + a1 * b1 + a2 + b2
    if d == 2:
        M = d2

    elif d == 3:
        d3 = b3 + a1 * b2 + a2 * b1 + a3
        M = d2 + d3

    elif d == 4:
        d3 = b3 + a1 * b2 + a2 * b1 + a3
        d4 = b4 + a1 * b3 + a2 * b2 + a3 * b1 + a4
        M = d2 + d3 + d4

    else:
        return decimal.Decimal(2) ** 300

    return 2 * M.log10() / log10Two + (3 * (k + 1 - f * mu)).log10() / log10Two - f * (
            1 - mu / beta).log10() / log10Two


# The cost of AGB over F_2
def AGBfor2(n, k, h):
    n = decimal.Decimal(n)
    k = decimal.Decimal(k)
    finalcost = 999999
    finalf = 99999
    finalmu = 999999
    for f in range(0, h + 1):
        for mu in range(0, math.floor(n / h) + 1):
            if f * mu < k + 1:
                subcost = subAGBfor2(n, k, h, f, mu)
                if finalcost > subcost:
                    finalcost = subcost
                    finalf = f
                    finalmu = mu

    # print(degree_conjfor2(n, k, h, finalf, finalmu))
    # print(finalf)
    # print(finalmu)
    # print("   for F2  n = " + str(n) + "     k = " + str(k) + "    h = " + str(h) + "  cost  " + str(round(finalcost, 8)))

    return round(finalcost, 2)


######################################################################### AGB 2.0 r2 & r3 #################################################


####################### The cost of AGB over F_q

# compute size of M1 and M2 for large field

def M1forq(n, k, h, fprime, muprime, degree, sizeM1):
    ff = fprime
    mumu = muprime
    n = decimal.Decimal(n)
    k = decimal.Decimal(k)
    h = decimal.Decimal(h)
    f = decimal.Decimal(fprime)
    mu = decimal.Decimal(muprime)
    beta = decimal.Decimal(math.floor(n / h))

    a1 = (beta - mu - 1) * f
    a2 = (beta - mu - 1) ** 2 * f * (f - 1) / 2
    a3 = (beta - mu - 1) ** 3 * f * (f - 1) * (f - 2) / 6
    a4 = (beta - mu - 1) ** 4 * f * (f - 1) * (f - 2) * (f - 3) / 24

    b1 = (beta - 1) * (h - f)
    b2 = (beta - 1) ** 2 * (h - f) * (h - f - 1) / 2
    b3 = (beta - 1) ** 3 * (h - f) * (h - f - 1) * (h - f - 2) / 6
    b4 = (beta - 1) ** 4 * (h - f) * (h - f - 1) * (h - f - 2) * (h - f - 3) / 24

    c1 = (h - 1)
    c2 = (h) * (h - 1) / 2
    c3 = (h + 1) * (h) * (h - 1) / 6
    c4 = (h + 2) * (h + 1) * (h) * (h - 1) / 24

    d2 = a1 + b1 + c1 + a1 * b1 + a1 * c1 + b1 * c1 + a2 + b2 + c2
    d3 = c3 + b1 * c2 + b2 * c1 + b3 + a1 * (b1 * c1 + b2 + c2) + a2 * (b1 + c1) + a3
    d4 = c4 + b1 * c3 + b2 * c2 + b3 * c1 + b4 + a1 * (b1 * c2 + b2 * c1 + b3 + c3) + a2 * (
            b2 + c2 + b1 * c1) + a3 * (b1 + c1) + a4

    sizeM1[ff][mumu][0] = d2
    sizeM1[ff][mumu][1] = d2 + d3
    sizeM1[ff][mumu][2] = d2 + d3 + d4

    if degree == 2:
        return sizeM1[ff][mumu][0]

    if degree == 3:
        return sizeM1[ff][mumu][1]

    if degree == 4:
        return sizeM1[ff][mumu][2]

    return decimal.Decimal(2) ** 300


def M2forq(n, k, h, fprime, muprime, degree, sizeM2):
    ff = fprime
    mumu = muprime
    n = decimal.Decimal(n)
    k = decimal.Decimal(k)
    h = decimal.Decimal(h)
    f = decimal.Decimal(fprime)
    mu = decimal.Decimal(muprime)
    beta = decimal.Decimal(math.floor(n / h))

    a1 = -(n - k - h)
    a2 = (n - k - h) * (n - k - h - 1) / 2
    a3 = -(n - k - h) * (n - k - h - 1) * (n - k - h - 2) / 6
    a4 = (n - k - h) * (n - k - h - 1) * (n - k - h - 2) * (n - k - h - 3) / 24

    b1 = (beta - mu - 1) * f
    b2 = (beta - mu - 1) ** 2 * f * (f - 1) / 2
    b3 = (beta - mu - 1) ** 3 * f * (f - 1) * (f - 2) / 6
    b4 = (beta - mu - 1) ** 4 * f * (f - 1) * (f - 2) * (f - 3) / 24

    c1 = (beta - 1) * (h - f)
    c2 = (beta - 1) ** 2 * (h - f) * (h - f - 1) / 2
    c3 = (beta - 1) ** 3 * (h - f) * (h - f - 1) * (h - f - 2) / 6
    c4 = (beta - 1) ** 4 * (h - f) * (h - f - 1) * (h - f - 2) * (h - f - 3) / 24

    d2 = a1 * b1 + a1 * c1 + b1 * c1 + a2 + b2 + c2
    d3 = c3 + b1 * c2 + b2 * c1 + b3 + a1 * (b1 * c1 + b2 + c2) + a2 * (b1 + c1) + a3
    d4 = c4 + b1 * c3 + b2 * c2 + b3 * c1 + b4 + a1 * (b1 * c2 + b2 * c1 + b3 + c3) + a2 * (b2 + c2 + b1 * c1) + a3 * (
            b1 + c1) + a4

    sizeM2[ff][mumu][0] = max(d2, Decimal(1))
    sizeM2[ff][mumu][1] = max(d2 + d3, Decimal(1))
    sizeM2[ff][mumu][2] = max(d2 + d3 + d4, Decimal(1))

    if degree == 2:
        return sizeM2[ff][mumu][0]

    if degree == 3:
        return sizeM2[ff][mumu][1]

    if degree == 4:
        return sizeM2[ff][mumu][2]

    return decimal.Decimal(2) ** 300


####################### The cost of AGB 2.0 r2 over F_q
# compute cost for parameter f and mu
def subOurAGB2forq(n, k, h, f, mu, sizeM1, sizeM2, finalcost):
    d = degree_conjforq(n, k, h, f, mu)
    beta = decimal.Decimal(math.floor(n / h))

    if d > 4:
        return 999999

    mincost = 999999
    minfprime = 99999
    minmuprime = 999999

    subgcost = 9999999999
    stepg = repeat * 3


    if sizeM2[f][mu][d - 2] > 0:
        M2 = decimal.Decimal(sizeM2[f][mu][d - 2])
    else:
        M2 = decimal.Decimal(M2forq(n, k, h, f, mu, d, sizeM2))

    #print("mincost  =  " + str(mincost)+"    M2  = "+str(M2))
    if M2 < 2:
        if sizeM1[f][mu][d - 2] > 0:
            M1 = decimal.Decimal(3 * (k + 1 - f * mu)) * decimal.Decimal(
                sizeM1[f][mu][d - 2]) ** Decimal(2)
        else:
            M1 = decimal.Decimal(3 * (k + 1 - f * mu)) * decimal.Decimal(
                M1forq(n, k, h, f, mu, d, sizeM1)) ** Decimal(2)

        subcost= M1 / ((1 - mu/ beta) ** f)
        mincost = subcost.log10()/ log10Two
        #print("mincost  =  "+str(mincost))

    for g in range(f * mu, -1, -1):

        if subgcost > mincost - decimal.Decimal(0.01):
            stepg = stepg - 1

        else:
            stepg = repeat * 3

        if stepg < 0:
            break

        mincost = min(mincost, subgcost)

        subgcost = 9999999999

        step = repeat

        fprime = 0
        muprime = 0
        for ff in range(f, 1, -1):


            if step < 0:
                break

            if g / ff > ff:
                break




            if g % ff == 0:
                fprime = ff
                muprime = math.floor(g / ff + 0.00001)




                if muprime > mu:
                    continue

                step = step - 1

                if fprime <= 1:
                    continue



                if sizeM1[fprime][muprime][d - 2] > 0:
                    M1 = decimal.Decimal(3 * (k + 1 - fprime * muprime)) * decimal.Decimal(
                        sizeM1[fprime][muprime][d - 2]) ** Decimal(2)
                else:
                    M1 = decimal.Decimal(3 * (k + 1 - fprime * muprime)) * decimal.Decimal(
                        M1forq(n, k, h, fprime, muprime, d, sizeM1)) ** Decimal(2)

                if sizeM2[fprime][muprime][d - 2] > 0:
                    M2 = decimal.Decimal(5.64) * decimal.Decimal(sizeM2[fprime][muprime][d - 2]) ** Decimal(2.8)
                else:
                    M2 = decimal.Decimal(5.64) * decimal.Decimal(
                        M2forq(n, k, h, fprime, muprime, d, sizeM2)) ** Decimal(2.8)


                mm = max(decimal.Decimal(sizeM1[fprime][muprime][d - 2]) - decimal.Decimal(sizeM2[fprime][muprime][d - 2]),1)

                M11 = decimal.Decimal(sizeM1[fprime][muprime][d - 2]) * mm ** Decimal(2) / mm.log10() * log10Two

                if M1 > M11:
                    M1 = M11

                if decimal.Decimal(sizeM2[fprime][muprime][d - 2]) > mm:
                    M2 = decimal.Decimal(5.64) * mm ** Decimal(2.8)


                M1 = max(1, M1)
                M2 = max(1, M2)
                T = M1 / ((1 - muprime / beta) ** fprime) + M2 / ((1 - mu / beta) ** f)


                subcost = T.log10() / log10Two

                subgcost = min(subgcost,subcost)
                if mincost > subcost:
                    minfprime = fprime
                    minmuprime = muprime

    if finalcost > mincost:
        '''
        with open("C:\\Users\\hans1\\Desktop\\AGB5.txt", "a") as text:
            text.write(str(datetime.datetime.now()) + "\n")
            text.write(
                "   for larger field   n = " + str(n) + "     k = " + str(
                    k) + "     d = " + str(
                    d) + "    h = " + str(h) + "    f = " + str(
                    f) + "    mu = "
                + str(mu) + "    minfprime = " + str(
                    minfprime) + "    minmuprime = "
                + str(minmuprime)  + "  finalcost  = " + str(round(finalcost, 8)) + "  mincost  " + str(
                    round(mincost, 8)) + "\n")
            text.write(
                " ==============================================7777====================================================  " + "\n")
            text.write("  " + "\n")
        '''

    return mincost


####################### The cost of AGB 2.0 r2 over F_q
# The cost of AGB over F_q
def OurAGB2forq(n, k, h):

    sizeM1 = zeors_array = np.zeros((h + 1, math.floor(n / h) + 1, 3))
    sizeM2 = zeors_array = np.zeros((h + 1, math.floor(n / h) + 1, 3))
    # print(size)

    searchg = search

    kk = k
    n = decimal.Decimal(n)
    k = decimal.Decimal(k)

    d0 = min(5, degree_conjforq(n, k, h, 0, 0))

    #print(d0)
    finalcost = 9999999

    if d0 == 5:
        finalcost = 9999999

    else:
        M1 = decimal.Decimal(3 * (k + 1)) * decimal.Decimal(
            M1forq(n, k, h, 0, 0, d0, sizeM1)) ** Decimal(2)
        finalcost = M1.log10() / log10Two

    if d0 == 2:
        '''
        with open("C:\\Users\\hans1\\Desktop\\AGB5.txt", "a") as text:
            text.write(str(datetime.datetime.now()) + "\n")
            text.write(
                "   for larger field   n = " + str(n) + "     k = " + str(
                    k) + "    h = " + str(
                    h) + "  f = 0   mu = 0   fprime=0 muprime = 0 fprimeprime=0 muprimeprime = 0  cost  " + str(
                    round(finalcost, 8)) + "\n")
            text.write(
                " ==============================================4444====================================================  " + "\n")
            text.write("  " + "\n")
        '''

        return round(finalcost, 2)

    finalf = 0
    finalmu = 0
    finalsubmu = 0

    subfcost = 999999
    submucost = 999999
    mincost = finalcost


    stepf = 1
    f = h+stepf

    endf = h - 4

    while 1 == 1:
        f = f - stepf

        subfcost = 999999
        # print("1")
        beginfinalmu = 0

        endfinalmu = math.floor(min((k + 1) / f, n / h - 1))

        #else:
        beginmu = 0
        endmu = math.floor(min(k / f, n / h - 1))

        if f == h:
            beginmu = -1
        else:
            beginmu = 0

        for mu in range(endmu, beginmu, -1):
            if degree_conjforq(n, k, h, f, mu) > 4:
                continue

            subcost = subOurAGB2forq(n, k, h, f, mu, sizeM1, sizeM2, mincost)



            if subcost < mincost:
                if subcost < submucost:
                    if subcost < subfcost:
                        finalf = f
                        finalmu = mu

            submucost = min(subcost, submucost)
            subfcost = min(submucost, subfcost)
            '''
            with open("C:\\Users\\hans1\\Desktop\\AGB5.txt", "a") as text:
                text.write(str(datetime.datetime.now()) + "\n")
                text.write(
                    "  ======1123123====== for larger field   n = " + str(n) + "     k = " + str(
                        k) + "    h = " + str(h) + "    endf = " + str(endf) + "    stepf = " + str(
                        stepf) + "    f = " + str(f) + "    mu = "
                    + str(mu) + "  cost  " + str(round(subcost, 8)) + "\n")
                text.write(
                    " ======================================9999===========================================================  " + "\n")
                text.write("  " + "\n")
            '''


        if stepf == 0:
            f = f - 1
            if f <= endf:
                break

        elif subfcost > mincost - decimal.Decimal(0.01):

            endf = f

            f = min( f + math.floor(stepf*1.5+0.00001),h)
            if stepf < 3:
                stepf = 0
            else:
                stepf = 1
        else:
            if stepf < 9:
                stepf = 2 * stepf

        mincost = min(mincost, subfcost)



    dd0 = degree_conjforq(n, k, h, finalf, finalmu)

    for f in range(endf,math.floor(endf/2),-1):
        mu = finalmu

        endmu = math.floor(min(k / f, n / h - 1))

        if mu > endmu:
            continue

        if degree_conjforq(n, k, h, f, mu) != dd0:
            break

        subcost = subOurAGB2forq(n, k, h, f, mu, sizeM1, sizeM2, mincost)
        searchg = searchg - 1
        '''
        with open("C:\\Users\\hans1\\Desktop\\AGB5.txt", "a") as text:
            text.write(str(datetime.datetime.now()) + "\n")
            text.write(
                "  ======   for larger field   n = " + str(n) + "     searchg = " + str(searchg)  + "     k = " + str(
                    k) + "    h = " + str(h) + "    f = " + str(f) + "    mu = "
                + str(mu) + "    dd0 = "
                + str(dd0) + "  cost  " + str(round(mincost, 8)) + "\n")
            text.write(
                " ======================================121299999============================================================  " + "\n")
            text.write("  " + "\n")
        '''
        if mincost > subcost:
            searchg = search
            endf = f
        mincost = min(mincost, subcost)


    for f in range(endf,0,-1):
        if searchg < 0:
            break

        for g in range(finalf * finalmu, kk+1, 1):
            if math.floor(g * h/n) >= f:
                break

            if g / f > f:
                break

            if g % f == 0:
                mu = math.floor(g / f + 0.00001)

                if mu > n/h - 1:
                    break

                searchg = searchg - 1

                if searchg < 0:
                    break

                subcost = subOurAGB2forq(n, k, h, f, mu, sizeM1, sizeM2, mincost)
                if mincost > subcost:
                    searchg = search

                '''
                if mincost > subcost:
                    with open("C:\\Users\\hans1\\Desktop\\AGB5.txt", "a") as text:
                        text.write(str(datetime.datetime.now()) + "\n")
                        text.write(
                            "  ======   for larger field   n = " + str(n) + "     searchg = " + str(searchg)  + "     k = " + str(
                                k) + "    h = " + str(h) + "    f = " + str(f) + "    mu = "
                            + str(mu) + "  mincost  " + str(round(mincost, 8)) + "\n")
                        text.write(
                            " ======================================99999============================================================  " + "\n")
                        text.write("  " + "\n")
                else:
                    with open("C:\\Users\\hans1\\Desktop\\AGB5.txt", "a") as text:
                        text.write(str(datetime.datetime.now()) + "\n")
                        text.write(
                            "=7777=====   for larger field  n = " + str(n) + "     searchg = " + str(searchg)  + "     k = " + str(
                                k) + "    h = " + str(h) + "    f = " + str(f) + "    mu = "
                            + str(mu) + "  mincost  " + str(round(mincost, 8)) + "\n")
                        text.write(
                            " ======================================99999============================================================  " + "\n")
                        text.write("  " + "\n")
                '''
                mincost = min(mincost, subcost)


    '''
    with open("C:\\Users\\hans1\\Desktop\\AGB5.txt", "a") as text:
        text.write(str(datetime.datetime.now()) + "\n")
        text.write(
            "   for larger field   n = " + str(n) + "     k = " + str(
                k) + "    h = " + str(h) + "    finalf = " + str(finalf) + "    finalmu = "
            + str(finalmu) + "  cost  " + str(round(mincost, 8)) + "\n")
        text.write(
            " ======================================3333============================================================  " + "\n")
        text.write("  " + "\n")
        text.write(
            " ======================================3333============================================================  " + "\n")
        text.write("  " + "\n")
        text.write(
            " ======================================3333============================================================  " + "\n")
        text.write("  " + "\n")
    '''

    return round(mincost, 2)


####################### The cost of AGB over F_2 =======================


# compute size of M1 and M2 for F2
def M1for2(n, k, h, fprime, muprime, degree, sizeM1):
    ff = fprime
    mumu = muprime

    n = decimal.Decimal(n)
    k = decimal.Decimal(k)
    h = decimal.Decimal(h)
    f = decimal.Decimal(fprime)
    mu = decimal.Decimal(muprime)
    beta = decimal.Decimal(math.floor(n / h))

    a1 = (beta - mu - 1) * f
    a2 = (beta - mu - 1) ** 2 * f * (f - 1) / 2
    a3 = (beta - mu - 1) ** 3 * f * (f - 1) * (f - 2) / 6
    a4 = (beta - mu - 1) ** 4 * f * (f - 1) * (f - 2) * (f - 3) / 24

    b1 = (beta - 1) * (h - f)
    b2 = (beta - 1) ** 2 * (h - f) * (h - f - 1) / 2
    b3 = (beta - 1) ** 3 * (h - f) * (h - f - 1) * (h - f - 2) / 6
    b4 = (beta - 1) ** 4 * (h - f) * (h - f - 1) * (h - f - 2) * (h - f - 3) / 24

    d2 = a1 + b1 + a1 * b1 + a2 + b2
    d3 = b3 + a1 * b2 + a2 * b1 + a3
    d4 = b4 + a1 * b3 + a2 * b2 + a3 * b1 + a4

    sizeM1[ff][mumu][0] = d2
    sizeM1[ff][mumu][1] = d2 + d3
    sizeM1[ff][mumu][2] = d2 + d3 + d4

    if degree == 2:
        return sizeM1[ff][mumu][0]

    if degree == 3:
        return sizeM1[ff][mumu][1]

    if degree == 4:
        return sizeM1[ff][mumu][2]

    return decimal.Decimal(2) ** 300


def M2for2(n, k, h, fprime, muprime, degree, sizeM2):
    ff = fprime
    mumu = muprime

    n = decimal.Decimal(n)
    k = decimal.Decimal(k)
    h = decimal.Decimal(h)
    f = decimal.Decimal(fprime)
    mu = decimal.Decimal(muprime)
    beta = decimal.Decimal(math.floor(n / h))

    a1 = (beta - mu - 1) * f
    a2 = (beta - mu - 1) ** 2 * f * (f - 1) / 2
    a3 = (beta - mu - 1) ** 3 * f * (f - 1) * (f - 2) / 6
    a4 = (beta - mu - 1) ** 4 * f * (f - 1) * (f - 2) * (f - 3) / 24

    b1 = (beta - 1) * (h - f)
    b2 = (beta - 1) ** 2 * (h - f) * (h - f - 1) / 2
    b3 = (beta - 1) ** 3 * (h - f) * (h - f - 1) * (h - f - 2) / 6
    b4 = (beta - 1) ** 4 * (h - f) * (h - f - 1) * (h - f - 2) * (h - f - 3) / 24

    c1 = -(n - k - 1)
    c2 = (n - k) * (n - k - 1) / 2
    c3 = -(n - k + 1) * (n - k) * (n - k - 1) / 6
    c4 = (n - k + 2) * (n - k + 1) * (n - k) * (n - k - 1) / 24

    d2 = a1 + b1 + c1 + a1 * b1 + a1 * c1 + b1 * c1 + a2 + b2 + c2
    d3 = c3 + b1 * c2 + b2 * c1 + b3 + a1 * (b1 * c1 + b2 + c2) + a2 * (b1 + c1) + a3
    d4 = c4 + b1 * c3 + b2 * c2 + b3 * c1 + b4 + a1 * (b1 * c2 + b2 * c1 + b3 + c3) + a2 * (b2 + c2 + b1 * c1) + a3 * (
            b1 + c1) + a4

    sizeM2[ff][mumu][0] = max(d2, Decimal(1))
    sizeM2[ff][mumu][1] = max(d2 + d3, Decimal(1))
    sizeM2[ff][mumu][2] = max(d2 + d3 + d4, Decimal(1))

    if degree == 2:
        return sizeM2[ff][mumu][0]

    if degree == 3:
        return sizeM2[ff][mumu][1]

    if degree == 4:
        return sizeM2[ff][mumu][2]

    return decimal.Decimal(2) ** 300


# compute cost for parameter f and mu

####################### The cost of AGB 2.0 r2 over F_2
# compute cost for parameter f and mu
def subOurAGB2for2(n, k, h, f, mu, sizeM1, sizeM2, finalcost):
    d = degree_conjfor2(n, k, h, f, mu)
    beta = decimal.Decimal(math.floor(n / h))

    if d > 4:
        return 999999

    mincost = 999999
    minfprime = 99999
    minmuprime = 999999

    subgcost = 9999999999
    stepg = repeat * 3
    for g in range(f * mu, 0, -1):

        if subgcost > mincost - decimal.Decimal(0.01):
            stepg = stepg -1

        else:
            stepg = repeat * 3

        if stepg < 0:
            break

        mincost = min(mincost,subgcost)

        subgcost = 9999999999

        step = repeat

        fprime = 0
        muprime = 0


        for ff in range(f, 0, -1):
            if step < 0:
                break

            if g / ff > ff:
                break

            if g % ff == 0:
                fprime = ff
                muprime = math.floor(g / ff + 0.00001)

                if muprime > mu:
                    continue

                step = step - 1

                if fprime <= 1:
                    continue

                if sizeM1[fprime][muprime][d - 2] > 0:
                    M1 = decimal.Decimal(3 * (k + 1 - fprime * muprime)) * decimal.Decimal(
                        sizeM1[fprime][muprime][d - 2]) ** Decimal(2)
                else:
                    M1 = decimal.Decimal(3 * (k + 1 - fprime * muprime)) * decimal.Decimal(
                        M1for2(n, k, h, fprime, muprime, d, sizeM1)) ** Decimal(2)

                if sizeM2[fprime][muprime][d - 2] > 0:
                    M2 = decimal.Decimal(5.64) * decimal.Decimal(sizeM2[fprime][muprime][d - 2]) ** Decimal(2.8)
                else:
                    M2 = decimal.Decimal(5.64) * decimal.Decimal(
                        M2for2(n, k, h, fprime, muprime, d, sizeM2)) ** Decimal(2.8)

                mm = max(decimal.Decimal(sizeM1[fprime][muprime][d - 2]) - decimal.Decimal(sizeM2[fprime][muprime][d - 2]),1)

                M11 = decimal.Decimal(sizeM1[fprime][muprime][d - 2]) * mm ** Decimal(2) / mm.log10() * log10Two

                if M1 > M11:
                    M1 = M11

                if decimal.Decimal(sizeM2[fprime][muprime][d - 2]) > mm:
                    M22 = decimal.Decimal(5.64) * mm ** Decimal(2.8)
                    M2 = M22

                M1 = max(1, M1)
                M2 = max(1, M2)

                T = M1 / ((1 - muprime / beta) ** fprime) + M2 / ((1 - mu / beta) ** f)
                subcost = T.log10() / log10Two

                subgcost = min(subgcost,subcost)
                if mincost > subcost:
                    minfprime = fprime
                    minmuprime = muprime


    if finalcost > mincost:
        '''
        with open("C:\\Users\\hans1\\Desktop\\AGB5.txt", "a") as text:
            text.write(str(datetime.datetime.now()) + "\n")
            text.write(
                "   for F2   n = " + str(n) + "     k = " + str(
                    k) + "    h = " + str(h) + "    f = " + str(
                    f) + "    mu = "
                + str(mu) + "    minfprime = " + str(
                    minfprime) + "    minmuprime = "
                + str(minmuprime) +  "  finalcost  = " + str(round(finalcost, 8)) + "  mincost  " + str(
                    round(mincost, 8)) + "\n")
            text.write(
                " ==============================================7777====================================================  " + "\n")
            text.write("  " + "\n")
        '''

    return mincost


####################### The cost of AGB 2.0 r2 over F_2
# The cost of AGB over F_2
def OurAGB2for2(n, k, h):
    sizeM1 = zeors_array = np.zeros((h + 1, math.floor(n / h) + 1, 3))
    sizeM2 = zeors_array = np.zeros((h + 1, math.floor(n / h) + 1, 3))
    # print(size)

    searchg = search

    kk = k
    n = decimal.Decimal(n)
    k = decimal.Decimal(k)

    d0 = min(5, degree_conjfor2(n, k, h, 0, 0))

    #print(d0)
    finalcost = 9999999

    if d0 == 5:
        finalcost = 9999999

    else:
        M1 = decimal.Decimal(3 * (k + 1)) * decimal.Decimal(
            M1for2(n, k, h, 0, 0, d0, sizeM1)) ** Decimal(2)
        finalcost = M1.log10() / log10Two


    if d0 == 2:
        '''
        with open("C:\\Users\\hans1\\Desktop\\AGB5.txt", "a") as text:
            text.write(str(datetime.datetime.now()) + "\n")
            text.write(
                "   for F2  n = " + str(n) + "     k = " + str(
                    k) + "    h = " + str(
                    h) + "  f = 0   mu = 0   fprime=0 muprime = 0 fprimeprime=0 muprimeprime = 0  cost  " + str(
                    round(finalcost, 8)) + "\n")
            text.write(
                " ==============================================4444====================================================  " + "\n")
            text.write("  " + "\n")
          '''


        return round(finalcost, 2)

    finalf = 0
    finalmu = 0
    finalsubmu = 0

    subfcost = 999999
    submucost = 999999
    mincost = finalcost


    stepf = 1
    f = h+stepf

    endf = h - 4

    while 1 == 1:
        f = f - stepf

        subfcost = 999999
        # print("1")
        beginfinalmu = 0

        endfinalmu = math.floor(min((k + 1) / f, n / h - 1))

        #else:
        beginmu = 0
        endmu = math.floor(min(k / f, n / h - 1))

        if f == h:
            beginmu = -1
        else:
            beginmu = 0


        for mu in range(endmu, beginmu, -1):
            if degree_conjfor2(n, k, h, f, mu) > 4:
                continue

            subcost = subOurAGB2for2(n, k, h, f, mu, sizeM1, sizeM2, mincost)



            if subcost < mincost:
                if subcost < submucost:
                    if subcost < subfcost:
                        finalf = f
                        finalmu = mu

            submucost = min(subcost, submucost)
            subfcost = min(submucost, subfcost)

            '''
            with open("C:\\Users\\hans1\\Desktop\\AGB5.txt", "a") as text:
                text.write(str(datetime.datetime.now()) + "\n")
                text.write(
                    "  ======1123123====== for F2   n = " + str(n) + "     k = " + str(
                        k) + "    h = " + str(h) + "    endf = " + str(endf) + "    stepf = " + str(
                        stepf) + "    f = " + str(f) + "    mu = "
                    + str(mu) + "  cost  " + str(round(subcost, 8)) + "\n")
                text.write(
                    " ======================================9999===========================================================  " + "\n")
                text.write("  " + "\n")
            '''


        if stepf == 0:
            f = f - 1
            if f <= endf:
                break

        elif subfcost > mincost - decimal.Decimal(0.01):

            endf = f

            f = min( f + math.floor(stepf*1.5+0.00001),h)
            if stepf < 3:
                stepf = 0
            else:
                stepf = 1
        else:
            if stepf < 9:
                stepf = 2 * stepf

        mincost = min(mincost, subfcost)



    dd0 = degree_conjfor2(n, k, h, finalf, finalmu)

    for f in range(endf,math.floor(endf/2),-1):
        mu = finalmu

        endmu = math.floor(min(k / f, n / h - 1))

        if mu > endmu:
            continue

        if degree_conjfor2(n, k, h, f, mu) != dd0:
            break

        subcost = subOurAGB2for2(n, k, h, f, mu, sizeM1, sizeM2, mincost)
        searchg = searchg - 1
        '''
        with open("C:\\Users\\hans1\\Desktop\\AGB5.txt", "a") as text:
            text.write(str(datetime.datetime.now()) + "\n")
            text.write(
                "  ======   for F2   n = " + str(n) + "     searchg = " + str(searchg)  + "     k = " + str(
                    k) + "    h = " + str(h) + "    f = " + str(f) + "    mu = "
                + str(mu) + "    dd0 = "
                + str(dd0) + "  cost  " + str(round(mincost, 8)) + "\n")
            text.write(
                " ======================================121299999============================================================  " + "\n")
            text.write("  " + "\n")
        '''
        if mincost > subcost:
            searchg = search
            endf = f
        mincost = min(mincost, subcost)


    for f in range(endf,0,-1):
        if searchg < 0:
            break

        for g in range(finalf * finalmu, kk+1, 1):
            if math.floor(g * h/n) >= f:
                break

            if g / f > f:
                break

            if g % f == 0:
                mu = math.floor(g / f + 0.00001)

                if mu > n/h - 1:
                    break

                searchg = searchg - 1

                if searchg < 0:
                    break

                subcost = subOurAGB2for2(n, k, h, f, mu, sizeM1, sizeM2, mincost)
                if mincost > subcost:
                    searchg = search
                '''
                if mincost > subcost:
                    with open("C:\\Users\\hans1\\Desktop\\AGB5.txt", "a") as text:
                        text.write(str(datetime.datetime.now()) + "\n")
                        text.write(
                            "  ======   for F2   n = " + str(n) + "     searchg = " + str(searchg)  + "     k = " + str(
                                k) + "    h = " + str(h) + "    f = " + str(f) + "    mu = "
                            + str(mu) + "  mincost  " + str(round(mincost, 8)) + "\n")
                        text.write(
                            " ======================================99999============================================================  " + "\n")
                        text.write("  " + "\n")
                else:
                    with open("C:\\Users\\hans1\\Desktop\\AGB5.txt", "a") as text:
                        text.write(str(datetime.datetime.now()) + "\n")
                        text.write(
                            "=7777=====   for F2 n = " + str(n) + "     searchg = " + str(searchg)  + "     k = " + str(
                                k) + "    h = " + str(h) + "    f = " + str(f) + "    mu = "
                            + str(mu) + "  mincost  " + str(round(mincost, 8)) + "\n")
                        text.write(
                            " ======================================99999============================================================  " + "\n")
                        text.write("  " + "\n")
                '''
                mincost = min(mincost, subcost)

    '''
    with open("C:\\Users\\hans1\\Desktop\\AGB5.txt", "a") as text:
        text.write(str(datetime.datetime.now()) + "\n")
        text.write(
            "   for F2   n = " + str(n) + "     k = " + str(
                k) + "    h = " + str(h) + "    finalf = " + str(finalf) + "    finalmu = "
            + str(finalmu) + "  cost  " + str(round(mincost, 8)) + "\n")
        text.write(
            " ======================================3333============================================================  " + "\n")
        text.write("  " + "\n")
        text.write(
            " ======================================3333============================================================  " + "\n")
        text.write("  " + "\n")
        text.write(
            " ======================================3333============================================================  " + "\n")
        text.write("  " + "\n")
    '''

    return round(mincost, 2)

##################### Table ###################################
def table():
    print("==================AGB2.0==================128==========larger=====================" + "\n")

    print(OurAGB2forq(1024, 652, 106))
    print(OurAGB2forq(1024 * 4, 1589, 172))
    print(OurAGB2forq(1024 * 16, 3482, 338))
    print(OurAGB2forq(1024 * 64, 7391, 667))
    print(OurAGB2forq(1024 * 256, 15336, 1312))
    print(OurAGB2forq(1024 * 1024, 32771, 2467))
    print(OurAGB2forq(1024 * 1024 * 4, 67440, 4788))

    print("==================AGB2.0==================128==========small=====================" + "\n")

    print(OurAGB2forq(1024, 652, 57))

    print(OurAGB2forq(1024 * 4, 1589, 98))
    print(OurAGB2forq(1024 * 16, 3482, 198))

    print(OurAGB2forq(1024 * 64, 7391, 389))
    print(OurAGB2forq(1024 * 256, 15336, 760))
    print(OurAGB2forq(1024 * 1024, 32771, 1419))
    print(OurAGB2forq(1024 * 1024 * 4, 67440, 2735))

    print("==================AGB2.0==================128==========larger=====================" + "\n")

    print(OurAGB2forq(1024, 652, 106))
    print(OurAGB2forq(1024 * 4, 1589, 172))
    print(OurAGB2forq(1024 * 16, 3482, 338))
    print(OurAGB2forq(1024 * 64, 7391, 667))
    print(OurAGB2forq(1024 * 256, 15336, 1312))
    print(OurAGB2forq(1024 * 1024, 32771, 2467))
    print(OurAGB2forq(1024 * 1024 * 4, 67440, 4788))

    print("==================AGB2.0==================2==========smaller=====================" + "\n")

    print(OurAGB2for2(1024, 652, 57))
    print(OurAGB2for2(1024 * 4, 1589, 98))
    print(OurAGB2for2(1024 * 16, 3482, 198))
    print(OurAGB2for2(1024 * 64, 7391, 389))
    print(OurAGB2for2(1024 * 256, 15336, 760))
    print(OurAGB2for2(1024 * 1024, 32771, 1419))
    print(OurAGB2for2(1024 * 1024 * 4, 67440, 2735))

    print("==================AGB2.0==================2==========larger=====================" + "\n")

    print(OurAGB2for2(1024, 652, 106))
    print(OurAGB2for2(1024 * 4, 1589, 172))
    print(OurAGB2for2(1024 * 16, 3482, 338))
    print(OurAGB2for2(1024 * 64, 7391, 667))

    print(OurAGB2for2(1024 * 256, 15336, 1312))
    print(OurAGB2for2(1024 * 1024, 32771, 2467))
    print(OurAGB2for2(1024 * 1024 * 4, 67440, 4788))

    print("==================AGB2.0======larger field==========Ferret==Wolverine============" + "\n")

    print(OurAGB2forq(642048, 19870, 2508))
    print(OurAGB2forq(10805248, 589760, 1319))

    print("==================AGB2.0======field F2=========Ferret==Wolverine============" + "\n")

    print(OurAGB2for2(609728, 36288, 1269))
    print(OurAGB2for2(10805248, 589760, 1319))




#####################      main() function   ###########################





def help():
    print(
            "============================================input error ================================================")
    print("input format of exact LPN: python  C:\\AGBscript.py n=1024 k=652 t=57  #(the concrete cost of AGB 2.0 to solve regular LPN over F2)")
    print("or python C:\\AGBscript.py n=1024 k=652 t=57 q #(the concrete cost of AGB 2.0 to solve regular LPN over larger field with field size |F|>2 )")
    print()

def main():
    if len(sys.argv) < 3 or len(sys.argv) > 6:
        help()
    elif 'n' in sys.argv[1] and 'k' in sys.argv[2] and 't' in sys.argv[3]:
        n = int(re.findall(r"\d+", sys.argv[1]).pop())
        k = int(re.findall(r"\d+", sys.argv[2]).pop())
        t = int(re.findall(r"\d+", sys.argv[3]).pop())

        if len(sys.argv) == 4:
            print("the concrete cost of AGB 2.0 to solve regular LPN (n=" + str(n) + ", k=" + str(k) + ", t=" + str(
                t) + ") over F2 :")
            print(OurAGB2for2(n, k, t))
            print()

        elif len(sys.argv) == 5 and 'q' in sys.argv[-1]:
            print("the concrete cost of AGB 2.0 to solve regular LPN (n=" + str(n) + ", k=" + str(k) + ", t=" + str(
                t) + ") over Fq for q>2 :")
            print(OurAGB2forq(n, k, t))

            print()

        else:
            help()

    else:
        help()

##################################
# Executed code
##################################
if __name__ == '__main__':
    table()
    main()
