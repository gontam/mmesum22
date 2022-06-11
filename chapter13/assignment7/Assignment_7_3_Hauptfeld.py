import matplotlib
import numpy as np
import matplotlib.pyplot as plt

def mc_pi(N):
    # Initialize variables for counting point numbers
    circle_count = 0
    square_count = 0
    pi = 0
    # Do the estimation N times
    for i in range(0, N):
        # Generate a random point
        x = np.random.uniform(-1, 1)
        y = np.random.uniform(-1, 1)
        # Calculate the distance between the point and (0,0)
        d = x**2 + y**2
        # Check if we're in the circle area
        if d <= 1:
            circle_count += 1
        # We're in the [-1, 1] square in all cases
        square_count += 1
        # Estimate pi using Monte Carlo
        approx_pi = 4 * circle_count / square_count
    return approx_pi

def mc_pi_stat(N, M):
    approx_pis = np.zeros(M)
    for i in range(0, M):
        approx_pis[i] = mc_pi(N)
    return np.mean(approx_pis), np.std(approx_pis, ddof=1)

def mc_pi_plt():
    # From 10^expo_from to 10^expo_to
    expo_from = 2
    expo_to = 5
    # Calculate total amount of iterations and initialize result arrays
    expo_count = expo_to - expo_from + 1
    sample_sizes = np.zeros(expo_count)
    mean_devs = np.zeros(expo_count)
    std_devs = np.zeros(expo_count)
    # Go through all the exponents
    for i in range(expo_count):
        print("Calculating N = 10^" + str(i+expo_from)+ "...")
        # Calculate for n^expo samples
        n = 10**(i+expo_from)
        mean_dev, std_dev = mc_pi_stat(n, 10)
        # Fill the results arrays
        sample_sizes[i] = np.log10(n)
        mean_devs[i] = mean_dev - np.pi
        std_devs[i] = np.log10(std_dev)
    # Plot the results
    fig, (ax1, ax2) = plt.subplots(1, 2)
    # Plot the mean deviation
    ax1.errorbar(sample_sizes, np.full(expo_count, np.pi), mean_devs)
    # Plot the standard deviation
    #plt.plot(sample_sizes, std_devs)
    p, V = np.polyfit(sample_sizes, std_devs, 1)
    ax2.plot(sample_sizes, std_devs, 'x')
    ax2.plot(sample_sizes, p * sample_sizes + V, '-')
    plt.show()

if __name__ == "__main__":
    mc_pi_plt()