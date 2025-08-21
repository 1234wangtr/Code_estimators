from lwyy.hardness_of_lpn import *
from hybrid.hybrid_quick import *

#####################      Bit security of LPN and dual LPN     ###########################
# We propose a non-asymptotic cost of the information set decoding algorithm, Pooled Gauss attack, and statistical decoding attack
# against the LPN problem over finite fields F_q or a ring Z_{2^\lambda} with parameters

# LPN parameters: N (number of queries),
#                k (length of secret),
#                t (Hamming weight of noise),
#                q (size of field) and
#                lambda (bit size of ring)

# dual LPN parametersï¼šn (corresponding to the number of COT/VOLE correlations),
#                     N (number of queries),
#                     t (Hamming weight of noise),
#                     q (size of field) and
#                     lambda (bit size of ring)

# noise parameters: exact or regular. Noise parameters can be either exact or regular. Default functions apply to exact noise, unless labeled as "regular".


def analysisfor2(N, k, t):
    # LWYY
    T1 = Gauss(N, k, t)
    print(f"Gauss={T1}")
    T2 = SDfor2(N, k, t)
    print(f"SD={T2}")
    T3 = SD2for2(N, k, t)
    print(f"SD2={T3}")
    T4 = SD_ISD(N, k, t)
    print(f"SDISD={T4}")
    T5 = BJMM_ISD(N, k, t)
    print(f"BJMM={T5}")

    # Hybrid

    return min(T1, T2, T3, T4, T5)


def analysisfor2regular(N, k, t):
    # LWYY
    T1 = AGBfor2(N, k, t)

    N1 = N - t
    k1 = k - t
    T2 = analysisfor2(N1, k1, t)

    # hybrid
    T3 = hybrid_2_quick(N,k,t)

    return min(T1, T2, T3)


def analysisfordual2(n, N, t):
    k = N - n
    return analysisfor2(N, k, t)


def analysisfordual2regular(n, N, t):
    k = N - n

    return analysisfor2regular(N, k, t)


def analysisforq(N, k, t, q):
    # LWYY
    T1 = SD_ISD_q(N, k, t, q)
    T2 = SDforq(N, k, t)
    T3 = Gauss(N, k, t)
    T4 = SD2forq(N, k, t, q)



    return min(T1, T2, T3, T4)


def analysisforqregular(N, k, t, q):
    # LWYY
    T1 = AGBforq(N, k, t)
    T2 = analysisforq(N, k, t, q)

    # hybrid
    T3 = hybrid_bigq_quick(N,k,t)

    return min(T1, T2, T3)


def analysisfordualq(n, N, t, q):
    k = N - n

    return analysisforq(N, k, t, q)


def analysisfordualqregular(n, N, t, q):
    k = N - n
    return analysisforqregular(N, k, t, q)


def analysisfor2lambda(N, k, t, lam):
    t = min(t, int(t * pow(2, lam - 1) / (pow(2, lam) - 1) + 1))

    return analysisfor2(N, k, t)


def analysisfor2lambdaregular(N, k, t, lam):
    return analysisfor2lambda(N, k, t, lam)


def analysisfordual2lambda(n, N, t, lam):
    k = N - n

    return analysisfor2lambda(N, k, t, lam)


def analysisfordual2lambdaregular(n, N, t, lam):
    k = N - n

    return analysisfordual2lambda(n, N, t, lam)