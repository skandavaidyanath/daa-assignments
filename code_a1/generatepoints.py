import random
import numpy as np


def generate(n):

    s = set()
    while len(s) < n:
        temp = (round(random.uniform(-1000, 1000), 50),
                round(random.uniform(-1000, 1000), 50))
        s.add(temp)

    with open("TestCases/points4.txt", "w+") as f:
        for x in s:
            f.write(str(x[0])+" "+str(x[1])+"\n")

    r = n//4
    x = np.arange(-r, r)
    y_1 = (r**2 - x**2)**.5
    y_2 = -(r**2 - x**2)**.5
    with open(f"TestCases/{n}_circle_points.txt", "w+") as f:
        for i in range(len(x)):
            f.write(str(x[i])+" "+str(y_1[i])+"\n")
            f.write(str(x[i])+" "+str(y_2[i])+"\n")


if __name__ == '__main__':
    generate(int(100))
