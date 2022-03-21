
def task1(input_string):
    input_string = input_string.lower()
    reversed_input_string = input_string[::-1]
    if input_string == reversed_input_string:
        print(reversed_input_string + ' True')
        return
    print(reversed_input_string + ' False')


inp_string = input('Input a string')
task1(inp_string)