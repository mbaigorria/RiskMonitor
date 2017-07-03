# -*- coding: utf-8 -*-

# Vencimientos

# Mayo y Diciembre (o Diciembre y Mayo) para maiz
# Mayo y Noviembre (o Noviembre y Mayo) para soja
# Mayo y Diciembre (o Diciembre (Z) y Mayo (K)) para trigo

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates

import numpy as np

import datetime
import calendar
import seaborn

import csv
import json

import os.path

import math
import time
import requests

def roundup(x):
	return int(math.ceil(x / 100.0)) * 100

def rounddown(x):
	return int(math.floor(x / 100.0)) * 100

def getMonthFromOptionCode(optionCode):
	codes = {'F': 'Jan', 'G': 'Feb', 'H': 'Mar', 'J': 'Apr', 'K': 'May', 
			'M': 'June', 'N': 'July', 'Q': 'Aug', 'U': 'Sept', 'V': 'Oct', 'X': 'Nov', 'Z': 'Dec'};
	if optionCode in codes:
		return codes[optionCode]
	else:
		return 'Invalid code.' 

def find_between(s, first, last):

	try:
		start = s.index( first ) + len( first )
		end = s.index( last, start )
		return s[start:end]
	except ValueError:
		print 'Invalid data.'
		exit()

def getPriceSeries(cropCode, monthCode, year, records=120):

	code = str(cropCode)+str(monthCode)+str(year)

	url = 'https://www.barchart.com/proxies/timeseries/queryeod.ashx?symbol='+code+'&data=daily&maxrecords='+str(records)
	headers = {'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_10_1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/39.0.2171.95 Safari/537.36'}

	response = requests.get(url, headers=headers)

	dates = []
	prices_open = []
	prices_high = []
	prices_low = []
	prices_close = []
	volumes = []

	for price_data in response.content.split('\n'):
		
		try:
			code, date, price_open, price_high, price_low, price_close, vol, _ = price_data.split(',')
		except:
			continue

		dates.append(date)
		prices_open.append(float(price_open))
		prices_high.append(float(price_high))
		prices_low.append(float(price_low))
		prices_close.append(float(price_close))
		volumes.append(float(vol))

	if len(dates) is 0:
		raise ValueError('failed to retrieve data')

	return dates, prices_close

def generateGraph(monthCode, assetType, expirationDate, underlyingPrice, thresholds):

	publishedMonths = {'soybean': ['X','K'], 'corn': ['K','Z'], 'wheat': ['K','Z']}

	if monthCode[0] not in publishedMonths[assetType]: return;

	day, month, year = expirationDate.split('-')

	if month == 'Dec':
		year = int(year) + 1
		year = str(year)

	cropCodes = {'soybean': 'ZS', 'wheat':'ZW', 'corn': 'ZC'}

	month = getMonthFromOptionCode(monthCode[0])

	date_list, price_list = getPriceSeries(cropCodes[assetType], monthCode[0], str(year)[-2:])

	monthToNumber = {v: k for k,v in enumerate(calendar.month_abbr)}

	day, month2, _ = expirationDate.split('-')
	expirationDate = year+'-'+str(monthToNumber[month2])+'-'+day

	assetTypeToSpanish = {'soybean': 'soja', 'wheat': 'trigo', 'corn': u'maíz'}

	# convert units to dollars per ton
	if assetType == 'corn':
		convertTometricTon = 39.3682;
	elif assetType == 'soybean':
		convertTometricTon = 36.7437;
	else:
		convertTometricTon = 36.7437;

	price_list = [p * convertTometricTon / 100.0 for p in price_list]
	thresholds = [p * convertTometricTon / 100.0 for p in thresholds]
	underlyingPrice = float(underlyingPrice) * convertTometricTon / 100.0

	# time series
	#sharey=True
	fig, axes = plt.subplots(ncols=2, gridspec_kw = {'width_ratios':[6, 1]})
	axes[0].plot_date(x=date_list, y=price_list, fmt="-", color="#FBBD16")
	axes[0].xaxis.set_major_formatter(mdates.DateFormatter('%d-%m-%Y'))
	axes[0].xaxis.set_major_locator(mdates.DayLocator(interval=len(date_list)/6))
	# axes[0].set_xlim(right=base+datetime.timedelta(days=10))
	axes[0].set_xlim(left=date_list[0], right=expirationDate)

	priceMargin = 20
	axes[0].set_ylim(top=max(max(price_list), max(thresholds)+priceMargin))
	axes[0].set_ylim(bottom=min(min(price_list), min(thresholds)-priceMargin))

	axes[0].set_ylabel(u'Dólares por tonelada de '+assetTypeToSpanish[assetType]+' al vencimiento')
	# axes[0].set_xlabel('Fecha')
	fig.autofmt_xdate()

	# stacked bar chart
	bar_width = 0.1
	y_max = max(price_list)

	colors = ['#7F8487', '#7BC6A2', '#63CCF0', '#7BC6A2','#7F8487']

	# thresholds = [0.2 for i in range(1,6)]
	percentages = ['10%','15%', '50%', '15%', '10%']

	# add stacked bars to graph
	bar_handles = []
	bottom = axes[0].get_ylim()[0]

	# stacked bar width
	width  = axes[1].get_xlim()[1]

	for i, d in enumerate(thresholds):
		bar_handles.append(axes[1].bar(0, d - bottom, width=width, color=colors[i%len(colors)], bottom=bottom))
		bottom = d

	bar_handles.append(axes[1].bar(0, axes[0].get_ylim()[1] - bottom, width=width, color=colors[(len(thresholds)+1)%len(colors)], bottom=bottom))

	# add bar labels
	for j in xrange(len(bar_handles)):
		for i, handle in enumerate(bar_handles[j].get_children()):
			bl = handle.get_xy()
			x = 0.55*handle.get_width() + bl[0]
			y = 0.40*handle.get_height() + bl[1]
			axes[1].text(x,y, "%s chance" % (percentages[j]), ha='center',fontsize=10)

	axes[1].set_xlim([0, bar_width])
	axes[1].xaxis.set_visible(False)
	axes[1].yaxis.tick_right()
	axes[1].set_ylim(axes[0].get_ylim())
	yticks = thresholds
	yticks = [min(x, axes[0].get_ylim()[1]) for x in thresholds]
	yticks.append(underlyingPrice)

	# yticks.append(axes[0].get_ylim()[0])
	axes[1].set_yticks(yticks, minor=False)


	plt.subplots_adjust(wspace=0, hspace=0)
	# plt.show()

	plt.savefig(assetType+'_'+month.lower()+'_'+year +'.png',bbox_inches='tight')

	print 'Graph generated! '+ assetType+'_'+month.lower()+'_'+year
	print thresholds

if __name__=="__main__":

	csvfile = open('../loadData.csv', 'r')
	csvreader = csv.reader(csvfile, delimiter=';')

	for row in csvreader:
		_, expirationDate, underlyingPrice, monthCode, assetType, finalJson, optionJson = row
		data = json.loads(finalJson)
		cut_points = data['indicator']['cut_points']
		thresholds = [price for price, density in cut_points]
		generateGraph(monthCode, assetType, expirationDate, underlyingPrice, thresholds)

	print 'Get any missing data from http://www.cmegroup.com/trading/agricultural/'