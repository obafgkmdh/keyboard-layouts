import sys
from collections import defaultdict
from itertools import combinations

alpha = "etoainshrlduymwcgf.pb,vk'jx;zq"

def sort3(i1, i2, i3):
    if i1 > i2: i1, i2 = i2, i1
    if i2 > i3: i2, i3 = i3, i2
    if i1 > i2: i1, i2 = i2, i1
    return alpha[i1] + alpha[i2] + alpha[i3]

if len(sys.argv) < 2:
    print(f"usage: {sys.argv[0]} input_file [corpus]")
    exit(0)

totalBigrams = 0
totalSkipgrams = 0
bigrams = defaultdict(int)
skipgrams = defaultdict(int)
if len(sys.argv) > 2:
    corpus = sys.argv[2]
    with open(corpus) as f:
        for line in f:
            prev = prev2 = 0
            for c in line:
                c = c.lower()
                if c == '"':
                    c = "'"
                elif c == ":":
                    c = ";"
                if c not in alpha:
                    prev2 = prev
                    prev = 0
                    continue
                if prev != 0:
                    totalBigrams += 1
                    bigrams[c, prev] += 1
                    bigrams[prev, c] += 1
                if prev2 != 0:
                    totalSkipgrams += 1
                    skipgrams[c, prev2] += 1
                    skipgrams[prev2, c] += 1
                prev2 = prev
                prev = c

with open(sys.argv[1]) as f:
    num = 0
    while 1:
        g1 = f.read(6)
        if not g1: break
        g2 = f.read(6)
        chars = set(alpha) - set(g1 + g2)
        g3 = []
        for i in range(5):
            g3.append(min(chars, key=alpha.index) + f.read(2))
            chars -= set(g3[-1])
        g3.append(sort3(*map(alpha.index, chars)))
        num += 1

        if totalBigrams > 0:
            sfb = sfs = 0
            for group in [g1, g2, *g3]:
                for i, j in combinations(group, 2):
                    sfb += bigrams[i, j]
                    sfs += skipgrams[i, j]
            print(f"{num:8d}:", g1, g2, *g3, f"sfb: {sfb/totalBigrams*100:.4f}%, sfs: {sfs/totalSkipgrams*100:.4f}%")
        else:
            print(f"{num:8d}:", g1, g2, *g3)
