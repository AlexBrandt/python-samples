
import numpy as np
import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt
import pylab as P
import sys

scores = []

for line in open(sys.argv[1]):
        scores.append(float(line.strip().split()[0]))

plt.hist(scores, bins=np.amax(scores)-np.amin(scores), color='blue')
plt.title("Histogram of Translated Protein Coding Sequence Lengths")
# plt.xlim([1,10])
plt.xlabel('Protein Length')
plt.ylabel('Counts')

txt = "Avg: %.2f\n" % (np.average(scores))
txt += "Std: %.2f\n" % (np.std(scores))
txt += "Min: %.2f\n" % (np.amin(scores))
txt += "Max: %.2f" % (np.amax(scores))

print txt

plt.text(plt.gca().get_xlim()[1] * .8, plt.gca().get_ylim()[1] * .8, txt)

plt.savefig(sys.argv[1].split('.')[0]+".pdf", bbox_inches='tight',alpha=.8)
