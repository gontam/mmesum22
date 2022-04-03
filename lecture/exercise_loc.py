class BinExercise:
    def __init__(self):
        self.inp_string = None
        ex_number = input('Which task shall be done? Use numerics from 1-4')
        if ex_number.isdigit():
            if ex_number == 1:
                inp_string = input('Input a string')
                self.task1(inp_string)
        else:
            print('Congrats, you failed!')
            exit()

    def task1(self, input_string):
        input_string = input_string.lower()
        reversed_input_string = input_string[::-1]
        if input_string == reversed_input_string:
            print(reversed_input_string + ' True')
            return
        print(reversed_input_string + ' False')

    @property
    def inp_string(self):
        return self.inp_string

    @inp_string.setter
    def inp_string(self, value):
        self.inp_string = value
