The files found in this folder are an attept to download the price series from the CME. The CME uses TradingView's chart widget to display these series.

To get the series, the client connects to the server. However, the data is not requested by HTTP in JSON format. The client opens a websocket to cmedata-idc.tradingview.com/socket.io/websocket, and through several messages and a protocol called MQTT, the data in JSON format is requested.

To start the connection, the correct origin host must be set (to replicate browser behaviour). This is not supported by most of the libraries I have seen, it's supported on websocket connections but not if the MQTT protocol is used on top of websockets.

To make this script work, my guess is you should first try to succesfully connect via websockets with a library that supports MQTT. There are libraries that support MQTT, but they don't really allow you to set the request host, meaning the target server will deny any connection. Then you should also be able to send normal websocket messages on the connection, without any MQTT package.

The relevant .js files are included in this folder. Feel free to play around with them!

To view the messages between server and client, just use any websocket sniffer on the CME website. Try to use Firefox's sniffer, since Chrome's default sniffer do not display MQTT packages.

http://www.cmegroup.com/trading/agricultural/grain-and-oilseed/soybean_quotes_globex_options.html

Let me know if you get this working!

Martin
