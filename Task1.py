
def isPalindrome(s):
    return s == s[::-1]

if __name__ == '__main__':

    # input
    string = str(input("Input: "))

    PalindromeFlag = isPalindrome(string)

    print(f"Output: {string} {isPalindrome(string)}")