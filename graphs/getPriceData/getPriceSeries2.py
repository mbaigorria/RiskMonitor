#!/usr/bin/env python3

import paho.mqtt.client as mqtt

# This is the Subscriber

def on_connect(client, userdata, flags, rc):
    print("Connected with result code "+str(rc))
    # client.subscribe("topic/motor-A/dt")

def on_message(client, userdata, msg):
	print(msg.payload)
    # if (msg.payload == 'Q'):
    #   m.stop()
    #   client.disconnect()
    # elif (-100 <= int(msg.payload) <= 100):
    #   m.duty_cycle_sp=msg.payload

client = mqtt.Client(transport="websockets")
client.connect('wss://cmedata-idc.tradingview.com/socket.io/websocket')

client.on_connect = on_connect
client.on_message = on_message

client.loop_forever()

