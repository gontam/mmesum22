#Write a script that separates a list into matrices


string = [[1, 3, 8, 9, 5, 1, 3, 7, 6], [9, 3, 4, 1, 5, 1, 6, 7, 2], [6, 1, 8, 7, 5, 4, 2, 9, 4], [7, 9, 8, 7, 5, 3, 2, 1, 4], [5, 1, 8, 7, 5, 3, 2, 1, 4], [4, 1, 3, 2, 9, 4, 8, 7, 5]]

#devide into n blocks
n = 3

print(f"length of list = {len(string)}")
print(f"length of first element of list = {len(string[1])}")
i=0
for x in string:
    #print(f"line {i}: {x}")
    i += 1

    if ((len(x)%n) != 0):
        print(f"CAUTION: length not divisible by {n}")
        continue
    else:
        # calc how loop var
        j = int(len(x) / n)
        for k in range(j):
            print(*(x[(k * n):((k + 1) * n)]))
