# Assignment 6
#Hein and Koranteng


import numpy as np
from skimage import io

DIM = 100
NGEN = int(input("Input number of generations to be simulated:  "))


grid = np.zeros((DIM,DIM), dtype=np.int8)

def pad_array(arr):
    var = np.copy(arr)
    var = np.pad(var,1, mode='constant')
    var = np.pad(var,1, mode='constant')
    return var


def unpad_array(arr):
    return arr[2:DIM+2,2:DIM+2]

def compare_array(arr1,arr2):
    cnt = 0
    for i in range(3):
        for j in range(3):
            if (arr1[i,j]==arr2[i,j]):
               cnt += 1

    if (cnt == 9):
        return True
    else: 
        return False

def simple_rule(nhood):
    lifes = np.sum(nhood)

    if(nhood[1,1]== 1):
        if(lifes < 3 or lifes > 5):
            return 0
        else: 
            return 1
    else: 
        if(lifes >= 3):
            return 1
        else: 
            return 0    

def apply_rule(padded_grid):
    for i in range(2,DIM):
        for j in range(2,DIM):   
            nhood = grid[i-2: i+1, j-2: j+1]
            padded_grid[i,j] = simple_rule(nhood)
    return padded_grid

def init_patten(grid):
    grid[1,1] = 1
    grid[1,2] = 1    
    grid[1,3] = 1
    grid[2,1] = 1

    grid[11,1] = 1
    grid[11,2] = 1
    grid[12,1] = 1
    grid[12,2] = 1

    grid[1,20] = 1
    grid[2,20] = 1
    grid[3,20] = 1
    grid[4,20] = 1
      
    grid[12,12] = 1
    grid[11,12] = 1    
    grid[12,13] = 1
    grid[12,11] = 1

    grid[20,20] = 1 
    grid[21,20] = 1
    grid[20,19] = 1 
    grid[19,21] = 1
    
 ###################   

    grid[61,61] = 1
    grid[61,62] = 1    
    grid[61,63] = 1
    grid[62,61] = 1

    grid[71,71] = 1
    grid[71,72] = 1
    grid[72,71] = 1
    grid[72,72] = 1

    grid[91,90] = 1
    grid[92,90] = 1
    grid[93,90] = 1
    grid[94,90] = 1
      
    grid[52,52] = 1
    grid[51,52] = 1    
    grid[52,53] = 1
    grid[52,51] = 1

    grid[10,70] = 1 
    grid[11,70] = 1
    grid[10,69] = 1 
    grid[ 9,71] = 1 

    return grid


#nhood = grid[:3,0:3]
#patt = np.array([[1,1,1],[0,0,0],[0,0,0]])

grid = init_patten(grid)
print(grid)

for i in range(NGEN):
    tmp = pad_array(np.copy(grid))
    im = (grid * 255).astype(np.uint8)
    fname = str(i)+ '.png'
    io.imsave(fname,im)
    grid = unpad_array(apply_rule(tmp))

    
            
print(grid)
