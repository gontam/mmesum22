import random
import numpy as np
import math
import matplotlib.pyplot as plt
#reference: https://www.geeksforgeeks.org/estimating-value-pi-using-monte-carlo/

#function for calculation estimated PI value
def mc_pi(N):
  circle_points = 0
  square_points = 0
  INTERVAL = N
  for i in range(INTERVAL ** 2):

    # Randomly generated x and y values from a
    # uniform distribution
    # Range of x and y values is -1 to 1
    rand_x = random.uniform(-1, 1)
    rand_y = random.uniform(-1, 1)

    # Distance between (x, y) from the origin
    origin_dist = rand_x ** 2 + rand_y ** 2

    # Checking if (x, y) lies inside the circle
    if origin_dist <= 1:
        circle_points += 1

    square_points += 1

    # Estimating value of pi, works correctly
    # pi= 4*(no. of points generated inside the
    # circle)/ (no. of points generated inside the square)
    pi = 4 * circle_points / square_points
  print(pi)
  return pi

def  mc_pi_stat(N,M):
    arr = []
    for i in range(M-1):
        arr.append(mc_pi(N))
#here we return the mean and standard deviation, works correctly
    print(np.mean(arr))
    print(np.std(arr,ddof=1))
    return np.mean(arr), np.std(arr,ddof=1)

def mc_pi_plt():
    n = 10
    n_arr = []
    m_arr = []
    std_arr = []
    tmp = []
#we generate the different sample sizes for the plot
    for i in range(5):
        n_arr.append(n*10)
        n = n*10
#we take the means and standard deviations
    for i in range(5):
        tmp.append(mc_pi_stat(n_arr[i], 10))

    print(tmp)
    m_arr,std_arr = tmp.pop()
    plt.figure()
    plt.errorbar(n_arr, m_arr, std_arr, linestyle='None', marker='^')
    plt.show()
#everything works well except tmp in def mc_pi_plt().

mc_pi_plt()