import sys
total = 0
for line in open(sys.argv[1]):
    total += int(line.strip())

print total
