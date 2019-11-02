import numpy as np
import seaborn as sns


def get_sigmoid():
    '''
    get sigmoid coordinates
    '''
    x = np.arange(-6,6,0.2)
    y = 1 / (1 + np.exp(-x))
    return x, y
    

def get_sin():
    '''
    get sin coordinates
    '''
    x = np.arange(-6,6,0.2)
    y = np.sin(x)
    return x, y


def make_file(x, y, filename):
    '''
    write points to a file
    '''
    with open(filename, 'w') as f:
        for (i, j) in zip(x, y):
            f.write(str(i) + ' ' + str(j) + '\n')


def plot(x, y):
    '''
    plots points given coordinate lists
    '''
    fig = sns.scatterplot(x, y).get_figure()
    fig.savefig('plot.png')


if __name__ == '__main__':
    x, y = get_sigmoid()
    make_file(x, y, 'data/sigmoid.txt')
    # plot(x, y)
    
    x, y = get_sin()
    make_file(x, y, 'data/sin.txt')
    # plot(x, y)
