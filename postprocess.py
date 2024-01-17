import sys

alpha = "etoainshrlduymwcgf.pb,vk'jx;zq"

def sort3(i1, i2, i3):
	if i1 > i2: i1, i2 = i2, i1
	if i2 > i3: i2, i3 = i3, i2
	if i1 > i2: i1, i2 = i2, i1
	return alpha[i1] + alpha[i2] + alpha[i3]

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
		print(f"{num:8d}:", g1, g2, *g3)