
word = input()

def is_palindrom(string):
    length = len(string)
    for i in range(0, length//2):
        if string[i] != string[length-1-i]:
            print(string, "False")
            return
    print(string, "True")
    return

is_palindrom(word)