def is_palindrome(value):
    return value[::-1] == value

test_string = input("input: ")
print(f"output: {test_string} {is_palindrome(test_string)}")
