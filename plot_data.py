#%%
#py_run_string(python_code)
import matplotlib.pyplot as plt

def plot_data(x, y):
    plt.plot(x, y)
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.title('Plot of X and Y Data')
    plt.show()

# Sample data
x = [1, 2, 3, 4, 5]
y = [2, 4, 6, 8, 10]


# Plotting data
plot_data(x, y)

# %%
