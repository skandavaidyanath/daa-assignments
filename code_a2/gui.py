import seaborn as sns
import tkinter as tk
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import subprocess


def get_data_points(filename):
    '''
    get two lists of x, y coordinates given a filename
    '''
    x, y = [], []
    with open(filename, 'r') as f:
        for line in f:
            currentX, currentY = line[:-1].split(' ')
            x.append(float(currentX))
            y.append(float(currentY))
    return x, y


def get_line_points(x, y, filename):
    '''
    get indices of points representing line endpoints
    '''
    x_line, y_line = [], []
    with open('out.txt', 'r') as f:
        for line in f:
            ind = int(line[:-1])
            x_line.append(x[ind])
            y_line.append(y[ind])
    return x_line, y_line


LARGE_FONT = ("Verdana", 12)


class main(tk.Tk):
    '''
    controller class
    '''

    def __init__(self, *args, **kwargs):

        tk.Tk.__init__(self, *args, **kwargs)

        tk.Tk.wm_title(self, "DAA A2 2018-19 SEM 2")

        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        page = Page(container, self)
        page.grid(row=0, column=0, sticky="nsew")
        self.show_frame(page)

    def show_frame(self, frame):
        frame.tkraise()


class Page(tk.Frame):
    '''
    class for the page
    '''

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        label = tk.Label(self, text='Segmented Linear Regression', font=LARGE_FONT)
        label.pack(pady=10, padx=10)

        slide = tk.Scale(self, from_=0, to=7.5, length=300,
                         resolution=0.01, orient=tk.HORIZONTAL)
        slide.pack()

        E = tk.Entry(self)
        E.pack()

        button = tk.Button(self, text="Plot",
                            command=lambda: self.plot(slide.get(), E.get()))
        button.pack()

        self.ax, self.fig = None, None
        self.canvas, self.toolbar = None, None


    def plot(self, c, filename):
        if self.ax is not None:
            self.ax.clear()
            self.fig.clf()
        
        # execute algorithm
        output = subprocess.check_output(['./a.out', f'{c}', f'data/{filename}.txt'])

        x_data, y_data = get_data_points(f'data/{filename}.txt')
        x_line, y_line = get_line_points(x_data, y_data, f'data/{filename}.txt')
        self.ax = sns.scatterplot(x_data, y_data)
        self.ax = sns.lineplot(x='x', y='y', data={'x': x_line, 'y': y_line}, ax=self.ax)
        self.fig = self.ax.get_figure()

        if self.canvas is None:
            self.canvas = FigureCanvasTkAgg(self.fig, self)
            self.canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

            self.toolbar = NavigationToolbar2Tk(self.canvas, self)
            self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.canvas.draw()
        self.toolbar.update()


app = main()
app.mainloop()
