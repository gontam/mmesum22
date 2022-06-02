# Implementation cellular automata with "Conway's Game of Life rule"
# Import necessary libraries:
import numpy as np
import argparse
import matplotlib.animation as animation
import matplotlib.pyplot as plt


# Randomized Grid for every cell:
def rndmGrd(m):
    return np.random.choice(vl, m * m, p=[.4, .6]).reshape(m, m)


# Glider for Cell(0,0)
def aglide(i, j, grd):
    glide = np.array([[0, 0, 1],
                     [1, 0, 1],
                     [0, 1, 1]])
    grd[i:i+3, j:j+3] = glide


# Glider Gun on Cell(0,0)
def aGlideG(i, j, grd):
    gn = np.zeros(11*38).reshape(11, 38)

    gn[5][1] = gn[5][2] = 1
    gn[6][1] = gn[6][2] = 1

    gn[3][13] = gn[3][14] = 1
    gn[4][12] = gn[4][16] = 1
    gn[5][11] = gn[5][17] = 1
    gn[6][11] = gn[6][15] = gn[6][17] = gn[6][18] = 1
    gn[7][11] = gn[7][17] = 1
    gn[8][12] = gn[8][16] = 1
    gn[9][13] = gn[9][14] = 1
    gn[1][25] = 1
    gn[2][23] = gn[2][25] = 1
    gn[3][21] = gn[3][22] = 1
    gn[4][21] = gn[4][22] = 1
    gn[5][21] = gn[5][22] = 1
    gn[6][23] = gn[6][25] = 1
    gn[7][25] = 1

    gn[3][35] = gn[3][36] = 1
    gn[4][35] = gn[4][36] = 1

    grd[i:i + 11, j:j + 38] = gn


def update(frame, img, grd, M):
    # Compute Neighbour sum:
    newGrd = grd.copy()
    for i in range(M):
        for j in range(M):
            ttl = int(grd[i, (j - 1) % M] + grd[i, (j + 1) % M] +
                      grd[(i - 1) % M, j] + grd[(i + 1) % M, j] +
                      grd[(i - 1) % M, (j - 1) % M] + grd[(i - 1) % M, (j + 1) % M] +
                      grd[(i + 1) % M, (j - 1) % M] + grd[(i + 1) % M, (j + 1) % M]) / 2

            # Conway rule:
            if grd[i, j] == 1:
                if(ttl < 2) or (ttl > 3):
                    newGrd[i, j] = 0
            else:
                if ttl == 3:
                    newGrd[i, j] = 1

    # Update
    img.set_data(newGrd)
    grd[:] = newGrd[:]
    return img


# Grid settings
on = 1
off = 0
vl = [on, off]

# Plot the Game of Life:
def plotting():
    pars = argparse.ArgumentParser(description="Welcome to Game of Life!")
    pars.add_argument('--grid-size', dest='N', required=False)
    pars.add_argument('--mov-file', dest='movfile', required=False)
    pars.add_argument('--interval', dest='interval', required=False)
    pars.add_argument('--glider', action='store_true', required=False)
    pars.add_argument('--gosper', action='store_true', required=False)
    args = pars.parse_args()

    m = int(input("Grid size: "))
    if args.N and int(args.N) > 8:
        m = int(args.N)

    u_int = int(input("Update Interval in seconds: "))
    if args.interval:
        u_int = int(args.interval)

    grd = np.array([])

    if args.glider:
        grd = np.zeros(m * m).reshape(m, m)
        aglide(1, 1, grd)
    elif args.gosper:
        grd = np.zeros(m * m).reshape(m, m)
        aGlideG(10, 10, grd)
    else:
        grd = rndmGrd(m)

    # Animation:
    fig, ax = plt.subplots()
    img = ax.imshow(grd, interpolation='nearest')
    anim = animation.FuncAnimation(fig, update, fargs=(img, grd, m), frames=10, interval=u_int, save_count=50)

    if args.movfile:
        anim.save(args.movfile, fps=60, extra_args=['-vcodec', 'libx264'])
    plt.show()


# main
def main():
    plotting()


# Main:
if __name__ == '__main__':
    main()
