import sys
import re

from all_estimator import *

def help():
    print(
        "============================================input error ================================================")
    print("input format of exact LPN: console.py N=1024 k=652 t=57 exact")
    print("============================================================================================")
    # print(
    #     "or console.py N=1024 k=652 t=57 lambda=13 exact  #(bit security of exact LPN with ring size 2^lambda)")
    print("or console.py N=1024 k=652 t=57 q=13 exact #(bit security of exact LPN with field size q")
    print("or console.py n=1024 N=4096 t=88 exact #(bit security of dual exact LPN)")
    # print(
    #     "or console.py n=1024 N=4096 t=88 lambda=13 exact #(bit security of dual exact LPN with ring size 2^lambda)")
    print("or console.py n=1024 N=4096 t=88 q=13 exact #(bit security of dual exact LPN with field size q")
    print(" ============================================================================================")
    print("input format of regular LPN: console.py N=1024 k=652 t=57 regular")
    # print(
    #     "or console.py N=1024 k=652 t=57 lambda=13 regular #(bit security of regular LPN with ring size 2^lambda)")
    print("or console.py N=1024 k=652 t=57 q=13 regular #(bit security of regular LPN with field size q")
    print("or console.py n=1024 N=4096 t=88 regular #(bit security of dual regular LPN)")
    # print(
    #     "or console.py n=1024 N=4096 t=88 lambda=13 regular #(bit security of dual regular LPN with ring size 2^lambda)")
    print("or console.py n=1024 N=4096 t=88 q=13 regular #(bit security of dual regular LPN with field size q")
    print(" ============================================================================================")
    # print("or console.py n=82 k=41 t=4 q=2 m=41")
    print()

#####################      main() function   ###########################

def main():
    # print(len(sys.argv))
    if len(sys.argv) < 4 or len(sys.argv) > 6:
        help()
    elif 'n' in sys.argv[1] and 'N' in sys.argv[2] and 't' in sys.argv[3]:
        n = int(re.findall(r"\d+", sys.argv[1]).pop())
        N = int(re.findall(r"\d+", sys.argv[2]).pop())
        t = int(re.findall(r"\d+", sys.argv[3]).pop())

        if len(sys.argv) == 5 and 'x' in sys.argv[-1]:
            print("bit security of dual exact LPN (n=" + str(n) + ", N=" + str(N) + ", t=" + str(t) + "):")
            print(analysisfordual2(n, N, t))
            print()

        elif len(sys.argv) == 5 and 'r' in sys.argv[-1]:
            print("bit security of dual regular LPN (n=" + str(n) + ", N=" + str(N) + ", t=" + str(t) + "):")
            print(analysisfordual2regular(n, N, t))
            print()

        elif 'q' in sys.argv[-2] and 'x' in sys.argv[-1]:
            q = int(re.findall(r"\d+", sys.argv[-2]).pop())
            print("bit security of dual exact LPN (n=" + str(n) + ", N=" + str(N) + ", t=" + str(t) + ", q=" + str(
                q) + "):")
            print(analysisfordualq(n, N, t, q))
            print()

        elif 'q' in sys.argv[-2] and 'r' in sys.argv[-1]:
            q = int(re.findall(r"\d+", sys.argv[-2]).pop())
            print(
                "bit security of regular LPN (n=" + str(n) + ", N=" + str(N) + ", t=" + str(t) + ", q=" + str(q) + "):")
            print(analysisfordualqregular(n, N, t, q))
            print()


        elif 'lambda' in sys.argv[-2] and 'x' in sys.argv[-1]:
            lam = int(re.findall(r"\d+", sys.argv[-2]).pop())
            print("bit security of dual exact LPN (n=" + str(n) + ", N=" + str(N) + ", t=" + str(t) + ", lambda=" + str(
                lam) + "):")
            print(analysisfordual2lambda(n, N, t, lam))
            print()

        elif 'lambda' in sys.argv[-2] and 'r' in sys.argv[-1]:
            lam = int(re.findall(r"\d+", sys.argv[-2]).pop())
            print(
                "bit security of dual regular LPN (n=" + str(n) + ", N=" + str(N) + ", t=" + str(t) + ", lambda=" + str(
                    lam) + "):")
            print(analysisfordual2lambdaregular(n, N, t, lam))
            print()

        else:
            help()


    elif 'N' in sys.argv[1] and 'k' in sys.argv[2] and 't' in sys.argv[3]:
        N = int(re.findall(r"\d+", sys.argv[1]).pop())
        k = int(re.findall(r"\d+", sys.argv[2]).pop())
        t = int(re.findall(r"\d+", sys.argv[3]).pop())

        if len(sys.argv) == 5 and 'x' in sys.argv[-1]:
            print("bit security of exact LPN (N=" + str(N) + ", k=" + str(k) + ", t=" + str(t) + "):")
            print(analysisfor2(N, k, t))
            print()

        elif len(sys.argv) == 5 and 'r' in sys.argv[-1]:
            print("bit security of regular LPN (N=" + str(N) + ", k=" + str(k) + ", t=" + str(t) + "):")
            print(analysisfor2regular(N, k, t))
            print()

        elif 'q' in sys.argv[-2] and 'x' in sys.argv[-1]:
            q = int(re.findall(r"\d+", sys.argv[-2]).pop())
            print("bit security of exact LPN (N=" + str(N) + ", k=" + str(k) + ", t=" + str(t) + ", q=" + str(
                q) + "):")
            print(analysisforq(N, k, t, q))
            print()


        elif 'q' in sys.argv[-2] and 'r' in sys.argv[-1]:
            q = int(re.findall(r"\d+", sys.argv[-2]).pop())
            print("bit security of regular LPN (N=" + str(N) + ", k=" + str(k) + ", t=" + str(t) + ", q=" + str(
                q) + "):")
            print(analysisforqregular(N, k, t, q))
            print()

        elif 'lambda' in sys.argv[-2] and 'x' in sys.argv[-1]:
            lam = int(re.findall(r"\d+", sys.argv[-2]).pop())
            print("bit security of exact LPN (N=" + str(N) + ", k=" + str(k) + ", t=" + str(t) + ", lambda=" + str(
                lam) + "):")
            print(analysisfor2lambda(N, k, t, lam))
            print()



        elif 'lambda' in sys.argv[-2] and 'r' in sys.argv[-1]:
            lam = int(re.findall(r"\d+", sys.argv[-2]).pop())
            print("bit security of regular LPN (N=" + str(N) + ", k=" + str(k) + ", t=" + str(t) + ", lambda=" + str(
                lam) + "):")
            print(analysisfor2lambdaregular(N, k, t, lam))
            print()

        else:
            help()


    else:
        help()

if __name__ == '__main__':
    main()