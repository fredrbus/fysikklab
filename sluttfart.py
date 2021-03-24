import math
import numpy as np

sluttfarter = {}

for n in range(1, 11):
    path = "./tracked_data/txy_" + str(n) + ".txt"
    f = open(path)
    lines = f.readlines()
    t1, x1, y1 = [float(e) for e in lines[-2].split("\t")]
    t2, x2, y2 = [float(e) for e in lines[-1].split("\t")]
    dt = t2 - t1
    ds = math.sqrt((x2 - x1)**2+(y2-y1)**2)
    v = ds / dt
    sluttfarter["forsøk " + str(n)] = np.round(v, 3)


verdier = list(sluttfarter.values())

sluttfarter_snitt = np.mean(verdier)

sluttfarter_sd = np.std(verdier)

sluttfarter_varians = np.var(verdier)

print("\nForsøk \t\t Sluttfart")
print("-"*27)
for n in range(1, 11):
    print(str(n) + " \t\t " + str(sluttfarter["forsøk " + str(n)]))


print("")
print("Snittfart:" + "\t" + str(sluttfarter_snitt))
print("Standardavvik:" + "\t" + str(sluttfarter_sd))
print("Varians:" + "\t" + str(sluttfarter_varians))

f.close()
