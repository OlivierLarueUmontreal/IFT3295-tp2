# verification de la Q2, pas à mettre dans le devoir!!!! # pseudogitignore

import os
S = 1
with open("similarite.csv") as f:
    for line in f:
        
        line = line.strip().split(",")
        summed = sum(list(map(lambda l: int(l),line)))
        print(f"S{S} score : {summed}")
        S += 1