import matplotlib.pyplot as plt
import numpy as np
import csv as csv

pi = np.pi

def readTSV(filename):
    dataset = []
    with open(filename) as fd:
        rd = csv.reader(fd, delimiter="\t")
        for row in rd:
            dataset.append(row)
    return dataset
        
s01Dataset = readTSV("s01.tsv")
s02Dataset = readTSV("s02.tsv")
s03Dataset = readTSV("s03.tsv")

def alpha1(S):
    return {
        0.1: 0.6981434804,
        0.2: 0.7114813170,
        0.3: 0.5548265538,
    }[S]
def alpha3(S):
    return {
        0.1: -0.3887653237,
        0.2: -0.1981685189,
        0.3: -0.03459089271,
    }[S]
def alpha5(S):
    return {
        0.1: 0.5733180778,
        0.2: 0.1968598368,
        0.3: 0.1911936493,
    }[S]
def beta1(S):
    return {
        0.1: -0.4062782223,
        0.2: -1.005090317,
        0.3: -0.9112266596,
    }[S]
def beta3(S):
    return {
        0.1: 0.3781478169,
        0.2: 1.257158136,
        0.3: 0.9981923712,
    }[S]
def beta5(S):
    return {
        0.1: -0.2027769333,
        0.2: -0.7457245580,
        0.3: -0.6894817537,
    }[S]

def c1(S, omega, V, mu, tau):
    return pi * (beta1(S) * omega + beta1(S) * omega * V**2 * mu * np.sin(omega * tau))

def c2(S, omega):
    return 0

def c3(S, omega):
    return 0

def solve(S):
    l = 0.6
    E = 2.1*10**11
    ro = 7800
    r1 = 0.01
    r2 = 0.0075
    s1 = 0.000314159
    s2 = 0.000176714
    m1 = 0.5
    m2 = 0.5
    # c1 = 10**6
    # c2 = 10**6
    h1 = 1
    h2 = 1
    J1 = 7.853981635*10**-9
    J2 = 2.485048877*10**-9
    mu = 0.6
    tau = 0.01

    alpha = 0.9975160950
    gamma = 8.959122501 * 10**5
    beta = 0.06/pi * np.sqrt(alpha * gamma)
    omega = np.sqrt(gamma/alpha)

    def v(x):
        return 1142.694438 * x**6 - 2056.849989 * x**5 + 1393.577770 * x**4 - 438.1833312 * x**3 + 51.60145250 * x**2 + 3.484168346 * x + 0.17221
    
    V = v(0.3)

    def c1():
        return pi * (beta * omega + beta1(S) * omega * V**2 + alpha1(S) * mu * np.sin(omega * tau))

    def c2():
        return 3/4 * pi * V**4 * (alpha3(S) * (1 - mu * np.cos(omega * tau)) ** 2 * mu * np.sin(omega * tau) + alpha3(S) * mu ** 3 * (np.sin(omega * tau)) ** 3 + beta3(S) * omega ** 3)

    def c3():
        return 5/8 * alpha5(S) * V**6 * (1 - mu * np.cos(omega * tau))**4 * mu * pi * np.sin(omega * tau) + 5/8 * alpha5(S) * V**6 * mu**5 * pi * np.sin(omega * tau) ** 5 + 5/4 * alpha5(S) * V**6 * (1 - mu * np.cos(omega * tau)) ** 2 * mu**3 * pi * np.sin(omega * tau) ** 3 + 5/8 * beta5(S) * V**6 * pi * omega**5 

    a1 = 0

    def a2():
        return 1/2 * (np.sqrt(-2 * c3() * (c2() - np.sqrt(c2() ** 2 - 4 * c3() * c1())))) / c3()

    def a3():
        return 1/2 * (np.sqrt(-2 * c3() * (c2() + np.sqrt(c2()**2 - 4 * c3() * c1())))) / c3()

    return [a1, a2(), a3()]

def solveAndPrintResFor(S):
    res = solve(S)
    print("Для s =", S)
    print("a1 =", res[0], ";\na2 =", res[1], ";\na3 =", res[2])

def solveAndPrintAllResults():
    solveAndPrintResFor(0.1)
    solveAndPrintResFor(0.2)
    solveAndPrintResFor(0.3)

def main():
    solveAndPrintAllResults()

if __name__ == "__main__":
    main()