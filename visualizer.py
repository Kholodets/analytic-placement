import numpy as np
import matplotlib.pyplot as plt

pref = open("pre_spreading.txt", "r")
postf = open("post_spreading.txt", "r")

prex = np.array([])
prey = np.array([])
for line in pref:
    elm = line.split()
    prex = np.append(prex, float(elm[1]) )
    prey = np.append(prey, float(elm[2]) )
    #print("x: ", float(elm[2]), "y: ", float(elm[4]))


postx = np.array([])
posty = np.array([])
for line in postf:
    elm = line.split()
    postx = np.append(postx, float(elm[1]))
    posty = np.append(posty, float(elm[2]))

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10,4))

ax[0].scatter(prex, prey, alpha = 0.25, s=1)
ax[0].set_title('pre-spreading')
ax[1].scatter(postx, posty, alpha = 0.25, s=1)
ax[1].set_title('post-spreading')

plt.savefig("visualizer_plot.png")


pref.close()
postf.close()
