def is_palindrome(value):
    return value[::-1] == value

test_string = input("Input: ")
print(f"Output: {is_palindrome(test_string)}")