import string

def is_palindrome(value):
    return value[::-1] == value

palindrome = str(input())
print(f"Output: {is_palindrome(palindrome)}")
