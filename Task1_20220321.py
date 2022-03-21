# Task 1 - palindrome check

# Write a function that checks whether a string is a palindrome
# Input must be provided by the user in the console
# Output must be presented in the console

def is_palindrome():
    # string input
    test_str = str(input("Input: "))

    # make it suitable for caseless comparison
    test_str = test_str.casefold()

    # reverse string
    reverse_str = reversed(test_str)

    # check if the string is equal to its reverse
    if list(test_str) == list(reverse_str):
        print("The string " + test_str+ " is a palindrome.")
    else:
        print("The string " + test_str+ " is not a palindrome.")

is_palindrome()