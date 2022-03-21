import math


def palindromCheck(textinput: str):
    palindrom = True;

    length = len(textinput);
    lengthHalf = math.floor(length/2);

    for i in range(lengthHalf-1):
        palindrom &= textinput[i]==textinput[length-1-i];

    return palindrom;