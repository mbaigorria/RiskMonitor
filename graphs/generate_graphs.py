import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates

import numpy as np
import datetime
import seaborn

numdays = 30
base = datetime.datetime.today()
date_list = [base - datetime.timedelta(days=x) for x in range(0, numdays)]
price_series = np.random.rand(30)

# time series
#sharey=True
fig, axes = plt.subplots(ncols=2, gridspec_kw = {'width_ratios':[6, 1]})
axes[0].plot_date(x=date_list, y=price_series, fmt="r-")
axes[0].xaxis.set_major_formatter(mdates.DateFormatter('%d-%m-%Y'))
axes[0].xaxis.set_major_locator(mdates.DayLocator(interval=5))
axes[0].set_xlim(right=base+datetime.timedelta(days=10))
axes[0].set_ylabel('Dolares por tonelada de soja al vencimiento')
axes[0].set_xlabel('Fecha')
fig.autofmt_xdate()

# stacked bar chart
bar_width = 0.1
y_max = max(price_series)

colors = 'rgcwm'


thresholds = [0.2 for i in range(1,6)]
percentages = ['10%','15%', '50%', '15%', '10%']

print thresholds
print percentages
print y_max

# add stacked bars to graph
bar_handles = []
bottom = 0
for i, d in enumerate(thresholds):
	bar_handles.append(axes[1].bar(0, d, width=bar_width, color=colors[i%len(colors)], bottom=bottom))
	bottom += d

# add bar labels
for j in xrange(len(bar_handles)):
	for i, handle in enumerate(bar_handles[j].get_children()):
		bl = handle.get_xy()
		x = 0.5*handle.get_width() + bl[0]
		y = 0.5*handle.get_height() + bl[1]
		axes[1].text(x,y, "%s chance" % (percentages[j]), ha='center',fontsize=10)


axes[1].set_xlim([0, bar_width])
axes[1].xaxis.set_visible(False)
axes[1].yaxis.tick_right()
yticks = [y_max * i * 0.2 for i in range(1, 6)]
yticks.append(0)
axes[1].set_yticks(yticks, minor=False)
# axes[1].set_ylabel('test')
# axes[1].yticks([y_max / i for i in range(1, 5)])
# axes[1].set_yticks([y_max / i for i in range(1, 5)], minor=False)
# axes[0].invert_xaxis()
plt.subplots_adjust(wspace=0, hspace=0)
# plt.show()


plt.savefig("sample.png",bbox_inches='tight')

# axes.set_xticks(ticks, minor=False)

# and

# axes.set_xticklabels(labels, fontdict=None, minor=False, **kwargs)

