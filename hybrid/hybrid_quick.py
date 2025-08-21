from math import *
import sys

refresh = True  # refresh:truncation  !refresh:without truncation


def bdt_hybrid(n, k, h):
    # bdt_hybrid_advanced(n, k, h)
    # print(f"***********")
    #hybrid_bigq(n, k, h)
    hybrid_bigq_quick(n,k,h)
    return hybrid_2_quick(n,k,h)
    #return bdt_hybrid_advanced(n, k, h)


def bdt_hybrid_test(half=True):
    '''
        The complexity of different parameters
    '''
    # half = True    # half:low weight parameters    !half:large weight parameters
    comp = []

    n, k, h = 2 ** 10, 652, 106
    if half:
        h = 57
    tot = bdt_hybrid(n, k, h)
    comp.append(tot)

    n, k, h = 2 ** 12, 1589, 172
    if half:
        h = 98
    tot = bdt_hybrid(n, k, h)
    comp.append(tot)

    n, k, h = 2 ** 14, 3482, 338
    if half:
        h = 198
    tot = bdt_hybrid(n, k, h)
    comp.append(tot)

    n, k, h = 2 ** 16, 7391, 667
    if half:
        h = 389
    tot = bdt_hybrid(n, k, h)
    comp.append(tot)

    n, k, h = 2 ** 18, 15336, 1312
    if half:
        h = 760
    tot = bdt_hybrid(n, k, h)
    comp.append(tot)

    n, k, h = 2 ** 20, 32771, 2467
    if half:
        h = 1419
    tot = bdt_hybrid(n, k, h)
    comp.append(tot)

    n, k, h = 2 ** 22, 64770, 4788
    if half:
        h = 2735
    tot = bdt_hybrid(n, k, h)
    comp.append(tot)

    comp.reverse()
    for cp in comp:
        print(cp)


def bdt_hybrid_advanced(n, k, h):
    beta = n // h
    # refresh parameter (with truncation)
    if refresh:
        r = n - k
        n = beta * h
        k = n - r

    print(f"n={n} k={k} h={h} n/h={n / h} beta={beta}")

    best = 1000
    best_u, best_f = -1, -1
    best_u0 = -1
    best_G1, best_G2, best_G3 = 0, 0, 0
    best_pr = 0
    best_trycomp = 0
    best_n0, best_m0 = 0, 0

    for u in range(beta - 1):
        for f in range(h):
            # \bar{f} each guess \bar{u}+1     h-\bar{f}-g each guess \bar{u}
            for q in range(h - f):  # represents the 'g' in paper
                # f u+1   h-f-q u
                if (n - h - f - (h - q) * u) <= 0:  # guess too many error
                    continue

                # debug

                log_pr = log2(1 - (u + 1) / beta) * f + log2(1 - (u) / beta) * (h - q - f)
                l = (n - k) // (beta - 1)

                n0 = k - h - (u + 1) * f - (h - f - q) * (u)
                u0 = 0
                G1 = (n - k - l * (beta - u - 2)) * (n - k) * ((n - k - l * (beta - u - 2)) + n0)
                G1 = max(G1, 1)
                log_G1 = log2(G1)
                log_G2 = 0
                log_G3 = 0
                try_comp = 0

                if n0 <= 0:  # no need to figure out quadratic equations
                    tot = -log_pr + log_G1
                else:
                    m0 = (h - f - q) * (beta - u - 2) * (beta - u - 1) / 2 + \
                         f * (beta - u - 2) * (beta - u - 3) / 2
                    w = (beta - 1) * q
                    u0 = w
                    v0 = n0 - w
                    if v0 < 0:
                        log_G3 = 1000
                        G2 = 2 ** 1000
                        G3 = 2 ** 1000
                    elif v0 == 0:
                        G2 = 1
                        G3 = 1
                        try_comp = log2(beta) * q
                    elif v0 <= sqrt(2 * m0) - 2:
                        G2 = v0 ** 4 * n0 ** 2 / 8
                        G2 /= log2(max(v0 ** 2 / 2, 2))
                        # G2 = (v0**2/2)**2.807
                        G2 = max(1, G2)

                        # G2 = 1

                        log_G2 = log2(G2)
                        G3 = 5 / 4 * (2 * v0 * v0 * w + v0 * w * w) + (2.5 * v0) ** 3 / log2(max(2.5 * v0, 2))
                        G3 = max(1, G3)
                        log_G3 = log2(G3)
                        try_comp = log2(beta) * q



                    else:
                        log_G3 = 1000
                        G2 = 2 ** 1000
                        G3 = 2 ** 1000
                    tot = -log_pr + log2(G1 + G2 + G3 * (2 ** try_comp))
                    # tot = -log_pr + max(log_G1, log_G2, log_G3 + try_comp)
                if tot < best:
                    best = tot
                    best_f = f
                    best_u = u
                    best_u0 = u0
                    best_n0 = n0
                    best_m0 = m0
                    best_G1 = log_G1
                    best_G2 = log_G2
                    best_G3 = log_G3
                    best_trycomp = try_comp
                    best_pr = log_pr

    f = best_f
    u = best_u
    u0 = best_u0
    q = u0 // (beta - 1)
    n0 = best_n0
    v0 = n0 - u0
    m0 = best_m0
    print(f"!!!best={best} !!!\n f={best_f} u={best_u} u0={best_u0} q={best_u0 // (beta - 1)}")
    print(f"m0={m0} n0={n0} v0={v0}")
    print(f"pr={best_pr} G1={best_G1} G2={best_G2} G3={best_G3} try_comp={best_trycomp}")
    # print concrete parameter

    print("----------")
    sys.stdout.flush()

    return best


debug = False


def hybrid_bigq(n, k, h):
    gauss_constant = 2.8

    beta = n // h
    # refresh parameter (with truncation)
    if refresh:
        r = n - k
        n = beta * h
        k = n - r

    print(f"n={n} k={k} h={h} n/h={n / h} beta={beta}")

    best = 1000
    best_u, best_f = -1, -1
    best_u0 = -1
    best_G1, best_G2, best_G3 = 0, 0, 0
    best_pr = 0
    best_trycomp = 0
    best_n0, best_m0 = 0, 0

    for u in range(0, beta - 1):
        for f in range(h):
            # \bar{f} each guess \bar{u}+1     h-\bar{f}-g each guess \bar{u}
            for q in range(h - f):  # represents the 'g' in paper
                # f u+1   h-f-q u
                if (n - f - (h - q) * u) <= 0:  # guess too many error
                    continue

                log_pr = log2(1 - (u + 1) / beta) * f + log2(1 - (u) / beta) * (h - q - f)
                l = (n - k) // (beta)

                n0 = k - (u + 1) * f - (h - f - q) * (u)
                m0 = (h - f - q) * (beta - u) * (beta - u - 1) / 2 + \
                     f * (beta - u - 2) * (beta - u - 1) / 2
                u0 = 0
                G1 = (n - k) ** 2 * n
                G1 = (n - k - l * (beta - u - 1)) * (n - k) * ((n - k - l * (beta - u - 1)) + n0)
                G1 = max(G1, 1)
                log_G1 = log2(G1)
                log_G2 = 0
                log_G3 = 0
                try_comp = 0

                if n0 <= 0:  # no need to figure out quadratic equations
                    tot = -log_pr + log_G1
                else:
                    m0 = (h - f - q) * (beta - u) * (beta - u - 1) / 2 + \
                         f * (beta - u - 2) * (beta - u - 1) / 2
                    w = (beta) * q
                    u0 = w
                    v0 = n0 - w
                    if v0 < 0:
                        log_G3 = 1000
                        G2 = 2 ** 1000
                        G3 = 2 ** 1000
                    elif m0 - v0 * (v0 + 1) // 2 >= 2.5 * (v0 + 2 * q + q * v0 + q * (q - 1) // 2):
                        remaining = v0 + 2 * q + q * v0 + q * (q - 1) // 2
                        G2 = (n0 * (n0 + 1) // 2 + n0) ** gauss_constant
                        # G2 /= log2(max((n0 * (n0 + 1) // 2 + n0), 2))
                        # G2 = (v0**2/2)**2.807
                        G2 = max(1, G2)
                        log_G2 = log2(G2)
                        # G3 = 5/4 * (2*v0*v0*w+w*w)+v0**3/log2(max(v0,2))
                        # G3 = max(1, G3, v0 * v0 * u0)
                        G3 = (2.5 * remaining) ** gauss_constant
                        # G3 /= log2(max(remaining, 2))
                        G3 += 2.5 * remaining * (v0 + 2 * w + w * v0 + w * (w - 1) // 2)
                        log_G3 = log2(G3)
                        try_comp = log2(beta) * q
                    else:
                        log_G3 = 1000
                        G2 = 2 ** 1000
                        G3 = 2 ** 1000
                    tot = -log_pr + log2(G1 + G2 + G3 * (2 ** try_comp))
                    # tot = -log_pr + max(log_G1, log_G2, log_G3 + try_comp)
                if tot < best:
                    best = tot
                    best_f = f
                    best_u = u
                    best_u0 = u0
                    best_n0 = n0
                    best_m0 = m0
                    best_G1 = log_G1
                    best_G2 = log_G2
                    best_G3 = log_G3
                    best_trycomp = try_comp
                    best_pr = log_pr

    f = best_f
    u = best_u
    u0 = best_u0
    q = u0 // (beta)
    n0 = best_n0
    v0 = n0 - u0
    m0 = best_m0
    remaining = v0 + q + q * v0 + q * (q - 1) // 2
    print(f"!!!best={best} !!!\n f={best_f} u={best_u} u0={best_u0} q={q}")
    print(f"m0={m0} n0={n0} v0={v0} remaining={remaining}")
    print(f"pr={best_pr} G1={best_G1} G2={best_G2} G3={best_G3} try_comp={best_trycomp}")
    # print concrete parameter

    print("----------")
    sys.stdout.flush()

    return best

def hybrid_bigq_quick(n, k, h, refresh=True):
    beta = n // h
    # refresh parameter (with truncation)
    if refresh:
        r = n - k
        n = beta * h
        k = n - r

    print(f"n={n} k={k} h={h} n/h={n / h} beta={beta}")

    best = float('inf')
    best_u, best_f = -1, -1

    # initial search range over u (integer)
    lo, hi = 0, beta - 1-1
    # evaluate endpoints
    f_lo = sub_bigq_quick_u(n, k, h, lo)
    f_hi = sub_bigq_quick_u(n, k, h, hi)
    evals = {lo: f_lo, hi: f_hi}

    # narrowing loop
    while hi - lo > 30:
        # split into 4 segments by 3 interior points
        d = (hi - lo) // 4
        m1 = lo + d
        m2 = lo + 2*d
        m3 = lo + 3*d

        for m in (m1, m2, m3):
            if m not in evals:
                evals[m] = sub_bigq_quick_u(n, k, h, m)

        #print(f"lo={lo}, m1={m1} f1={evals[m1]}, m2={m2} f2={evals[m2]}, m3={m3} f3={evals[m3]}, hi={hi}")

        # collect points and values
        us = [lo, m1, m2, m3, hi]
        fs = [evals[u] for u in us]
        # find index of minimal f
        idx = min(range(len(fs)), key=lambda i: fs[i])

        # shrink interval to neighbors of minimal index
        if idx == 0:
            # minimum at lo -> keep [lo, m1]
            hi = m1
            f_hi = evals[m1]
            evals[hi] = f_hi
        elif idx == len(us)-1:
            # at hi -> keep [m3, hi]
            lo = m3
            f_lo = evals[m3]
            evals[lo] = f_lo
        else:
            # at some mid -> keep [us[idx-1], us[idx+1]]
            lo = us[idx-1]
            hi = us[idx+1]
            f_lo = evals[lo]
            f_hi = evals[hi]
            evals[lo] = f_lo
            evals[hi] = f_hi

    # brute-force small window
    start = max(lo - 10, 0)
    end = min(hi + 10, beta - 1)
    for u in range(start, min(end + 1,beta-1-1)):
        val = evals.get(u, sub_bigq_quick_u(n, k, h, u))
        if val < best:
            best = val
            best_u = u

    print(f"T={best} u={best_u}")
    return best

def sub_bigq_quick_u(n,k,h,u):
    gauss_constant = 2.8
    beta = n // h
    best = 1000
    start = 0
    end = h
    Tstart = sub_bigq_quick_u_q(n, k, h, u, start)
    Tend = sub_bigq_quick_u_q(n, k, h, u, end)

    tmp_min = Tstart

    while end - start > 30:
        mid1 = int((end - start) / 3) + start
        mid2 = end - int((end - start) / 3)

        Tmid1 = sub_bigq_quick_u_q(n, k, h, u, mid1)
        Tmid2 = sub_bigq_quick_u_q(n, k, h, u, mid2)

        if Tmid1 > Tmid2:
            start = mid1
            tmp_min = Tmid2
        else:
            end = mid2
            tmp_min = Tmid1

    if start < 10:
        start = 0
    else:
        start = start - 10

    for q in range(start, min(end + 10, h), 1):
        T = sub_bigq_quick_u_q(n, k, h, u,q)
        # if u==6:
        #     print(f"u={u} f={f} T={T}")
        if T <= tmp_min:
            tmp_min = T

    best = min(tmp_min, best)
    return best

def sub_bigq_quick_u_q(n,k,h,u,q):

    gauss_constant = 2.8
    beta = n // h
    best = 1000
    for f in range(h - q):  # represents the 'g' in paper
        # f u+1   h-f-q u
        if (n - f - (h - q) * u) <= 0:  # guess too many error
            continue

        log_pr = log2(1 - (u + 1) / beta) * f + log2(1 - (u) / beta) * (h - q - f)
        l = (n - k) // (beta)


        n0 = k - (u + 1) * f - (h - f - q) * (u)
        m0 = (h - f - q) * (beta - u) * (beta - u - 1) / 2 + \
             f * (beta - u - 2) * (beta - u - 1) / 2
        u0 = 0
        G1 = (n - k) ** 2 * n
        G1 = (n - k - l * (beta - u - 1)) * (n - k) * ((n - k - l * (beta - u - 1)) + n0)
        G1 = max(G1, 1)
        log_G1 = log2(G1)
        log_G2 = 0
        log_G3 = 0
        try_comp = 0

        if n0 <= 0:  # no need to figure out quadratic equations
            tot = -log_pr + log_G1
        else:
            m0 = (h - f - q) * (beta - u) * (beta - u - 1) / 2 + \
                 f * (beta - u - 2) * (beta - u - 1) / 2
            w = (beta) * q
            u0 = w
            v0 = n0 - w
            if v0 < 0:
                log_G3 = 1000
                G2 = 2 ** 1000
                G3 = 2 ** 1000
            elif m0 - v0 * (v0 + 1) // 2 >= 2.5 * (v0 + 2 * q + q * v0 + q * (q - 1) // 2):
                remaining = v0 + 2 * q + q * v0 + q * (q - 1) // 2
                G2 = (n0 * (n0 + 1) // 2 + n0) ** gauss_constant

                G2 = max(1, G2)
                G3 = (2.5 * remaining) ** gauss_constant
                G3 += 2.5 * remaining * (v0 + 2 * w + w * v0 + w * (w - 1) // 2)
                try_comp = log2(beta) * q
            else:
                log_G3 = 1000
                G2 = 2 ** 1000
                G3 = 2 ** 1000
            tot = -log_pr + log2(G1 + G2 + G3 * (2 ** try_comp))
        if tot < best:
            best = tot

    return best

def hybrid_2_quick(n, k, h,refresh=True):
    beta = n // h
    # refresh parameter (with truncation)
    if refresh:
        r = n - k
        n = beta * h
        k = n - r

    print(f"n={n} k={k} h={h} n/h={n / h} beta={beta}")

    best = float('inf')
    best_u, best_f = -1, -1

    # initial search range over u (integer)
    lo, hi = 0, beta - 1-1
    # evaluate endpoints
    f_lo = sub_2_quick_u(n, k, h, lo)
    f_hi = sub_2_quick_u(n, k, h, hi)
    evals = {lo: f_lo, hi: f_hi}

    # narrowing loop
    while hi - lo > 30:
        # split into 4 segments by 3 interior points
        d = (hi - lo) // 4
        m1 = lo + d
        m2 = lo + 2*d
        m3 = lo + 3*d

        for m in (m1, m2, m3):
            if m not in evals:
                evals[m] = sub_2_quick_u(n, k, h, m)

        # print(f"lo={lo}, m1={m1} f1={evals[m1]}, m2={m2} f2={evals[m2]}, m3={m3} f3={evals[m3]}, hi={hi}")

        # collect points and values
        us = [lo, m1, m2, m3, hi]
        fs = [evals[u] for u in us]
        # find index of minimal f
        idx = min(range(len(fs)), key=lambda i: fs[i])

        # shrink interval to neighbors of minimal index
        if idx == 0:
            # minimum at lo -> keep [lo, m1]
            hi = m1
            f_hi = evals[m1]
            evals[hi] = f_hi
        elif idx == len(us)-1:
            # at hi -> keep [m3, hi]
            lo = m3
            f_lo = evals[m3]
            evals[lo] = f_lo
        else:
            # at some mid -> keep [us[idx-1], us[idx+1]]
            lo = us[idx-1]
            hi = us[idx+1]
            f_lo = evals[lo]
            f_hi = evals[hi]
            evals[lo] = f_lo
            evals[hi] = f_hi

    # brute-force small window
    start = max(lo - 10, 0)
    end = min(hi + 10, beta - 1)
    for u in range(start, min(end + 1,beta-1-1)):
        val = evals.get(u, sub_2_quick_u(n, k, h, u))
        if val < best:
            best = val
            best_u = u

    print(f"T={best} u={best_u}")
    return best

def sub_2_quick_u(n,k,h,u):
    beta = n // h
    best = 1000
    start = 0
    end = h
    Tstart = sub_2_quick_u_q(n, k, h, u, start)
    Tend = sub_2_quick_u_q(n, k, h, u, end)

    tmp_min = Tstart

    while end - start > 30:
        mid1 = int((end - start) / 3) + start
        mid2 = end - int((end - start) / 3)

        Tmid1 = sub_2_quick_u_q(n, k, h, u, mid1)
        Tmid2 = sub_2_quick_u_q(n, k, h, u, mid2)

        if Tmid1 > Tmid2:
            start = mid1
            tmp_min = Tmid2
        else:
            end = mid2
            tmp_min = Tmid1

    if start < 10:
        start = 0
    else:
        start = start - 10

    for q in range(start, min(end + 10, h), 1):
        T = sub_2_quick_u_q(n, k, h, u,q)
        if T <= tmp_min:
            tmp_min = T

    best = min(tmp_min, best)
    return best

def sub_2_quick_u_q(n,k,h,u,q):
    beta = n // h
    best = 1000
    for f in range(h - q):  # represents the 'g' in paper
        # f u+1   h-f-q u
        if (n - h - f - (h - q) * u) <= 0:  # guess too many error
            continue

        log_pr = log2(1 - (u + 1) / beta) * f + log2(1 - (u) / beta) * (h - q - f)
        l = (n - k) // (beta - 1)

        n0 = k - h - (u + 1) * f - (h - f - q) * (u)
        u0 = 0
        G1 = (n - k - l * (beta - u - 2)) * (n - k) * ((n - k - l * (beta - u - 2)) + n0)
        G1 = max(G1, 1)
        log_G1 = log2(G1)
        log_G2 = 0
        log_G3 = 0
        try_comp = 0

        if n0 <= 0:  # no need to figure out quadratic equations
            tot = -log_pr + log_G1
        else:
            m0 = (h - f - q) * (beta - u - 2) * (beta - u - 1) / 2 + \
                 f * (beta - u - 2) * (beta - u - 3) / 2
            w = (beta - 1) * q
            u0 = w
            v0 = n0 - w
            if v0 < 0:
                log_G3 = 1000
                G2 = 2 ** 1000
                G3 = 2 ** 1000
            elif v0 == 0:
                G2 = 1
                G3 = 1
                try_comp = log2(beta) * q
            elif v0 <= sqrt(2 * m0) - 2:
                G2 = v0 ** 4 * n0 ** 2 / 8
                G2 /= log2(max(v0 ** 2 / 2, 2))
                # G2 = (v0**2/2)**2.807
                G2 = max(1, G2)

                # G2 = 1

                log_G2 = log2(G2)
                G3 = 5 / 4 * (2 * v0 * v0 * w + v0 * w * w) + (2.5 * v0) ** 3 / log2(max(2.5 * v0, 2))
                G3 = max(1, G3)
                log_G3 = log2(G3)
                try_comp = log2(beta) * q



            else:
                log_G3 = 1000
                G2 = 2 ** 1000
                G3 = 2 ** 1000
            tot = -log_pr + log2(G1 + G2 + G3 * (2 ** try_comp))
            # tot = -log_pr + max(log_G1, log_G2, log_G3 + try_comp)
        if tot < best:
            best = tot
    return best

if __name__ == '__main__':

    # hybrid_2_quick(609728, 36288, 1269)
    # hybrid_2_quick(10805248, 589760, 1319)

    #hybrid_bigq_quick(642048, 19870, 2508)
    #hybrid_bigq_quick(10805248, 589760, 1319)

    bdt_hybrid_test(half=False)
    bdt_hybrid_test(half=True)
    #
    # hybrid_bigq(642048, 19870, 2508)
    # hybrid_bigq(10805248, 589760, 1319)
    #
    # bdt_hybrid_advanced(609728, 36288, 1269)
    # bdt_hybrid_advanced(10805248, 589760, 1319)

