axwidth = 3
axlength = 12
fontsize=32

plt.rc('axes',linewidth=axwidth)  #make sure to do this before making any figures

fig = plt.figure(figsize=(13,13))
ax = plt.gca()
<do plotting>

xtickloc = [list of xticks]
xtickstr = ['$'+str(kk)+'$' for kk in xtickloc]
ytickloc = [list of yticks]
ytickstr = ['$'+str(kk)+'$' for kk in ytickloc]

plt.xticks(xtickloc,xtickstr)
plt.yticks(ytickloc,ytickstr)

for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
for line in ax.get_xticklines() + ax.get_yticklines():
    line.set_markersize(18)
    line.set_markeredgewidth(3)
for tick in ax.xaxis.get_minor_ticks():
    tick.label1.set_fontsize(fontsize/2)
for tick in ax.yaxis.get_minor_ticks():
    tick.label1.set_fontsize(fontsize/2)

ax.tick_params(which='major',width=axwidth,length=axlength+5)
ax.tick_params(which='minor',width=axwidth,length=axlength)

plt.savefig()