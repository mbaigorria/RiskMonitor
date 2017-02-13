import asyncio
import websockets
import random
import logging

logger = logging.getLogger('websockets')
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler())

def randomHash():
	hash_string = ""
	for i in range(12):
		digit = round(60*random.random())
		hash_string = hash_string + "0123456789abcdefghijklmopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"[digit]
	return hash_string

def _sendRequest(websocket, a, b):
	c = {'m': a, 'p': b}
	d = str(c).replace(" ", "").replace("'", "\"")
	payload = "~m~" + str(len(d)) + "~m~" + d
	print("> " + payload)
	return websocket.send(payload)

async def hello():

	origin  = "https://s.tradingview.com"

	# seems these headers are not really required, use extra_headrs to add them via the connect method.
	headers = {'Host': 'cmedata-idc.tradingview.com',
				'User-Agent': 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:51.0) Gecko/20100101 Firefox/51.0',
				'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
				'Accept-Language': 'en-US,en;q=0.5',
				'Accept-Encoding': 'gzip, deflate, br',
				'Connection': 'keep-alive, Upgrade',
				'Pragma': 'no-cache',
				'Cache-Control': 'no-cache',
				'Upgrade:': 'websocket'}

	async with websockets.connect('wss://cmedata-idc.tradingview.com/socket.io/websocket', origin=origin) as websocket:

		chartSession  = "cs_" + randomHash()
		quoteSession  = "qs_" + randomHash()
		quoteSession2 = "qs_" + randomHash()

		response = await websocket.recv()
		print("< {}".format(response))

		_sendRequest(websocket, "chart_create_session", [chartSession,""])
		_sendRequest(websocket, "quote_create_session", [quoteSession])
		_sendRequest(websocket, "quote_set_fields", [quoteSession,"ch","chp","lp","fractional","minmov","minmove2","original_name","pricescale","pro_name","update_mode","volume","symbol_status","description","exchange","short_name","current_session","type","is_tradable"])
		_sendRequest(websocket, "quote_create_session", [quoteSession2])
		_sendRequest(websocket, "switch_timezone", [chartSession, "exchange"])
		_sendRequest(websocket, "resolve_symbol", [chartSession, "symbol_1","CBOT_GBX:ZSH2017"])
		_sendRequest(websocket, "create_series", [chartSession, "s1","s1","symbol_1","D",300])

		response = await websocket.recv()
		print("< {}".format(response))
		response = await websocket.recv()
		print("< {}".format(response))
		response = await websocket.recv()
		print("< {}".format(response))
		response = await websocket.recv()
		print("< {}".format(response))

		_sendRequest(websocket, "quote_add_symbols", [quoteSession,"CBOT:ZSH2017",{"flags":["force_permission"]}])
		_sendRequest(websocket, "create_study", [chartSession,"st1","st1","s1","Volume@tv-basicstudies-49",20,False])
		_sendRequest(websocket, "create_study", [chartSession,"st2","st1","s1","Sessions@tv-basicstudies-49"])

		response = await websocket.recv()
		print("< {}".format(response))

asyncio.get_event_loop().run_until_complete(hello())