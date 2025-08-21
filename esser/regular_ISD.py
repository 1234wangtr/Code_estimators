from math import log2, ceil, floor
from scipy.special import binom as binomial
import math


def representation_based_concrete_cost_depth_2(n, k, w):
    """
    Concrete cost of Representation-based Regular ISD (Section A.1) with depth 2
    Rounding issues are not considered
    The function outputs a list with two elements: the first element corresponds to the
    optimal parameters for the algorithm, the second element gives the minimum costs (both time and space complexities)
    The function minimizes only the time complexity, the space complexity is derived accordingly
    """

    beta = n // w
    r = n - k
    n = beta * w
    k = n - r
    b = n // w

    k_prime = k - w  # add parity-checks
    min_cost = [1000000000000000000000000000000, 0]
    best_params = [0, 0, 0, 0, 0, 0]

    # the algorithm raises an error in case the chosen ranges are not valid
    p_max = 40
    eps_x_max = 32

    for p in range(0, p_max, 4):
        for eps_x in range(0, eps_x_max, 2):
            p_x = p / 2 + eps_x  # weight

            v_approx = k_prime / w
            L1_approx = binomial(round(w / 2), p_x / 2) * k_prime ** (p_x / 2)
            if L1_approx <= 0:
                continue
            if w / 2 < p / 2:
                continue

            R = (binomial(p / 2, p / 4) * binomial(w / 2 - p / 2, eps_x / 2)) ** 2
            ell_approx = ceil(2 * log2(L1_approx))
            ell_min = floor(ell_approx * 0.5)
            ell_max = floor(ell_approx * 1.5)

            for ell in range(ell_min, ell_max):
                v = (k_prime + ell) / w

                # Num repr
                if R <= 0:
                    break
                ell_x = floor(log2(R))
                L1 = log2(binomial(round(w / 2), p_x / 2)) + log2(v) * (p_x / 2)  # list size, first level

                L_x1 = L1 * 2 - ell_x
                N_x = L_x1 * 2 - (ell - ell_x)

                # success probability
                if w / 2 < p / 2 or v >= b:
                    continue
                p_iter = (log2(binomial(floor(w / 2), round(p / 2))) +
                          log2(binomial(ceil(w / 2), round(p / 2))) +
                          log2((v / b)) * p +
                          log2(1 - v / b) * (w - p))

                # cost of one iteration
                T_iter = log2(n) + max(log2(n - k_prime) * 2, 2 + L1, 1 + L_x1, N_x)

                # overall cost
                cost = T_iter - p_iter

                # update min cost if new cost is smaller
                if cost < min_cost[0]:
                    min_cost = [cost, max(L1, L_x1) + log2(n)]
                    best_params = {"p": p, "eps_x": eps_x, "ell_x": ell_x,
                                   "ell_min": ell_min, "ell": ell, "ell_max": ell_max}

    if best_params["p"] == p_max:
        print("!!!!ERROR: p_max is too low for repr-based")

    if best_params["eps_x"] == eps_x_max:
        print("!!!!ERROR: eps_x_max is too low for repr-based")

    return min_cost

def permutation_based_concrete_cost_bigq(n, k, w):
    """
    Concrete cost of permutation-based Regular ISD (Section A.1)
    Rounding issues are not considered
    The function outputs a list with two elements, time complexity and space complexity
    """
    beta = n // w
    r = n - k
    n = beta * w
    k = n - r
    b = n // w

    k_prime = k #- w  # add parity-checks

    # success probability
    p_iter = (1 - k_prime/n)**w

    # cost of one iteration
    T_iter = (n - k_prime)**2 * n

    # overall cost
    # cost = log2(T_iter) - log2(p_iter)

    cost = log2(T_iter) - w*log2(1-k_prime/n)

    return [cost, log2(n - k) + log2(n)]

def permutation_based_concrete_cost(n,k,w):
    """
    Concrete cost of permutation-based Regular ISD (Section A.1)
    Rounding issues are not considered
    The function outputs a list with two elements, time complexity and space complexity
    """
    beta = n // w
    r = n - k
    n = beta * w
    k = n - r
    b = n // w

    k_prime  = k-w; #add parity-checks

    #success pr
    p_iter = (1-k_prime/n)**w;

    #cost of one iteration
    T_iter = (n-k_prime)**2*n;

    #overall cost
    cost = log2(T_iter)-log2(p_iter);

    return [cost, log2(n-k)+log2(n)];

def enumeration_based_concrete_cost(n, k, w):
    """
    Concrete cost of enumeration-based Regular ISD (Section A.1)
    Rounding issues are not considered
    The function outputs a list with two elements: the first element corresponds to the
    optimal parameters for the algorithm, the second element gives the minimum costs (both time and space complexities)
    The function minimizes only the time complexity, the space complexity is derived accordingly
    """

    beta = n // w
    r = n - k
    n = beta * w
    k = n - r
    b = n // w

    k_prime = k - w  # add parity-checks

    min_cost = [1000000000000000000000000000000, 0]
    best_params = []

    p_max = 30  # increase if p_max is found to be too low
    for p in range(0, p_max, 2):

        # we first find a suitable range for the optimal value of ell
        L1_approx = max(
            1,
            math.comb(round(w / 2), p // 2) * (k_prime / w) ** (p / 2)
        )
        ell = math.log2(L1_approx)
        ell_min = math.ceil(ell * 0.5)
        ell_max = min(round(ell * 1.5), n - k_prime)

        for ell in range(ell_min, ell_max):

            v = (k_prime + ell) / w  # number of coordinates per block

            # success probability
            if w / 2 < p / 2 or v / b >= 1:
                continue

            p_iter = (
                    math.log2(math.comb(math.floor(w / 2), round(p / 2)))
                    + math.log2(math.comb(math.ceil(w / 2), round(p / 2)))
                    + math.log2(v / b) * p
                    + math.log2(1 - v / b) * (w - p)
            )

            # cost of one iteration
            L = (
                    math.log2(math.comb(round(w / 2), p // 2))
                    + math.log2(v) * (p / 2)
            )
            T_iter = (
                    math.log2(n)
                    + max(math.log2(n - k_prime) * 2, 1 + L, L * 2 - ell)
            )

            # overall cost
            cost = T_iter - p_iter

            # update min cost  if new cost is smaller
            if cost < min_cost[0]:
                min_cost = [
                    cost,
                    max(L, math.log2(n - k)) + math.log2(n)
                ]  # time and space complexities
                best_params = {"p": p, "ell": ell}

    if best_params["p"] == p_max:
        print("!!!!ERROR: p_max is too low for enum-based")

    return best_params, min_cost


if __name__ == '__main__':
    n,k,w = 2**14,3482,338
    h = w
    beta = n // h
    refresh = True
    # refresh parameter (with truncation)
    if refresh:
        r = n - k
        n = beta * h
        k = n - r
    b = n//w
    print(representation_based_concrete_cost_depth_2(n,k,w))
    print(enumeration_based_concrete_cost(n, k, w))

    print(permutation_based_concrete_cost(n, k, w))