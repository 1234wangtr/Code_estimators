from math import *

def optimized_GRS(n,k,w,q,m):
    res = (n-k)**3
    res *= m**3
    res *= q**(w*ceil((k+1)*m/n) - m )
    res = log2(res)
    return res

def GRS(n,k,w,q,m):
    if n >= m:
        res = (n - k) ** 3
        res *= m ** 3
        res *= q ** (w * ceil(k * m / n) )
        res = log2(res)
    else:
        res = (n - k) ** 3
        res *= m ** 3
        res *= q ** (w * k)
        res = log2(res)
    return res

# def grobner_for_rank(n,k,w,q,m):
#     omega = 2.8
#     r = w
#     N = 0
#     i = 1
#     while i<=k:
#         N+= comb(n-i,r) * comb(k+b-1-i,b-1)
#         i+=1

def num_test1():
    n,k,w,q,m = 82,41,4,2,41
    res1 = GRS(n,k,w,q,m)
    res2 = optimized_GRS(n,k,w,q,m)
    print(f"res1={res1} res2={res2}")

    n, k, w, q, m = 106, 53, 5, 2, 53
    res1 = GRS(n,k,w,q,m)
    res2 = optimized_GRS(n,k,w,q,m)
    print(f"res1={res1} res2={res2}")

    n, k, w, q, m = 50, 32, 3, 2, 50
    res1 = GRS(n,k,w,q,m)
    res2 = optimized_GRS(n,k,w,q,m)
    print(f"res1={res1} res2={res2}")

    n, k, w, q, m = 112, 80, 4, 2, 112
    res1 = GRS(n,k,w,q,m)
    res2 = optimized_GRS(n,k,w,q,m)
    print(f"res1={res1} res2={res2}")


if __name__ == '__main__':
    num_test1()