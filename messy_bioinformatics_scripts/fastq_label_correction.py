import sys, gzip

f1 = gzip.open(sys.argv[1], 'rb')
f2 = gzip.open(sys.argv[2], 'wb')

broken_lines = [1, 3]

counter = 1

for line in f1:
    if counter % 4 in broken_lines:
        corrected_line = line[0:len(line)-2] + '2\n'
        f2.write(corrected_line)
    else:
        f2.write(line)
    counter += 1
