from FUSS import polplot as pp
import matplotlib.pyplot as plt

fig = plt.figure(1, figsize=(13, 6))
g=[]
loc=[131, 132, 133]
vel=[20000, 18000, 15000]
for i in range(len(loc)):
    grid = pp.axis(fig, loc = loc[i], phot_vel=vel[i], ang_grid='l', rad_grid = 'ul', vel_lim=[0,40000])
    g.append(grid)

h = pp.data(15,30, 35000)
he = pp.data(20,40,22000)
ca = pp.data(120,130,33000)

g[1].plot(h[0], h[1], lw=3, label='H I')
g[1].plot(he[0], he[1], lw=3, label='He I')
g[1].plot(ca[0], ca[1], lw=3, label='Ca II')

g[1].legend(bbox_to_anchor=(-3.5, 1.035, 5., .102), ncol=3)
plt.savefig('test_plot')
plt.show()

