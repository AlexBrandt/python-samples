import sys
sum=0
for line in open(sys.argv[1]):
    sum += int(line.strip().split()[0])
print sum
