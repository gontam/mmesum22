## Line Vector for storing data in it:
lines = []
## Reading in the whole data and store it in the vector "lines":
with open('G:\\My Drive\\Universität\\Master\\Master - 2. Semester (SS22) - Auslandssemester\\Bioinformatics\\BD137219.1.fasta','r') as f:
    lines = f.readlines()

## Getting information how long the fasta file is:
lengthFile = len(lines)
lastLine = len(lines[lengthFile - 1])

## If the last line only consists of two chars, put them in one before.
if(lastLine == 2):
    print("Hi")
else:
    print("no")

x = 1
stringNew = []
while x < lastLine:
    string = lines[x]
    string = string.strip('\n')
    stringNew += string
    print('String',x,'Verarbeitet.')
    print(string)
    x += 1


## Laboratory
string = lines[lengthFile - 2]
string = string.strip('\n')
stringNew = string + lines[lengthFile - 1]
print(string)
print(stringNew)
