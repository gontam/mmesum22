#!/usr/bin/env python3

# define a lambda function to check if string is a palindrome
is_palindrom = lambda string: string + " " + "True" if string == string[::-1] else string + " " + "False"

def main():
    # apply the defined function to a user input after asking the user for some string
    is_palindrom(input("Please enter some string: "))
    

if __name__ == '__main__':
    main()

