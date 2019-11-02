# importing the required module
import matplotlib.pyplot as plt

# x axis values
x1 = [8, 1000, 100, 10000, 10000]
x2 = [4, 13, 16, 29, 10000]

# corresponding y axis values
y1 = [0, 0, 0, 0.15, 0.15]

y2 = [0, 0, 0, 0.15, 3.86]

y3 = [0, 0, 0, 0.09, 0.34]

# plotting the points
plt.plot(x2, y1, label='Graham')
plt.plot(x2, y2, label='Jarvis')
plt.plot(x2, y3, label='Kirk Patrik Seidel')

# naming the x axis
# plt.xlabel('No of hull points')
plt.xlabel('no of hull points')
# naming the y axis
plt.ylabel('Elapsed seconds')

plt.legend()
# giving a title to my graph
plt.title('Time complexity')

# function to show the plot
fig1 = plt.gcf()
fig1.set_size_inches(18.5, 10.5)
plt.show()
fig1.savefig("time vs no of hull points.png", dpi=100)
