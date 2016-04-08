import sys

file = open(sys.argv[1])

for i in range(2):
    current = file.readline()

while current != "":
    print current,
    for i in range(4):
        current = file.readline()
