import math

numberOfPhiBins = 72
phiMin = 0
phiMax = 6.30
# numberOfPhiBins = 8
# phiMin = 0
# phiMax = 0.7

step = (phiMax - phiMin) / numberOfPhiBins

iPhi = list(range(0, numberOfPhiBins))
phi =[((v + 0.5) * step) for v in iPhi]

sin_phi = [round(math.sin(p), 5) for p in phi]
cos_phi = [round(math.cos(p), 5) for p in phi]

print "sin_phi is", sin_phi
print "cos_phi is", cos_phi
