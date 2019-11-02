import matplotlib.pyplot as plt


def plot(All_points, Hull_points):
    with open("bounds.txt") as f:
        content = f.readline()
    line = list(map(float, content.split()))
    plt.axis([line[0], line[2], line[1], line[3]])

    with open(All_points) as f:
        content = f.readlines()
    x = []
    y = []
    for l in content:
        line = list(map(float, l.split()))
        if(len(line) >= 2):
            x.append(line[0])
            y.append(line[1])
    plt.plot(x, y, 'ro')

    with open(Hull_points) as f:
        content = f.readlines()
    x = []
    y = []
    for l in content:
        line = list(map(float, l.split()))
        x.append(line[0])
        y.append(line[1])
    x.append(x[0])
    y.append(y[0])
    plt.plot(x, y, 'xb-')
    fig1 = plt.gcf()
    fig1.set_size_inches(18.5, 10.5)
    plt.show()
    fig1.savefig("k.png", dpi=100)


if __name__ == '__main__':
    plot("TestCases/points4.txt", "TestCases/outKS.txt")
