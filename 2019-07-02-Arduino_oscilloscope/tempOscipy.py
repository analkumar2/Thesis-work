# """
# exec(open('Arduino/tempOscipy.py').read())
# ldr.py
# Display analog data from Arduino using Python (matplotlib)
# Author: Mahesh Venkitachalam
# Website: electronut.in
# """
#
# import sys, serial, argparse
# import numpy as np
# from time import sleep
# from collections import deque
#
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation
#
#
# # plot class
# class AnalogPlot:
#   # constr
#   def __init__(self, strPort, maxLen):
#       # open serial port
#       self.ser = serial.Serial(strPort, 9600)
#
#       self.ax = deque([0.0]*maxLen)
#       self.ay = deque([0.0]*maxLen)
#       self.maxLen = maxLen
#
#   # add to buffer
#   def addToBuf(self, buf, val):
#       if len(buf) < self.maxLen:
#           buf.append(val)
#       else:
#           buf.pop()
#           buf.appendleft(val)
#
#   # add data
#   def add(self, data):
#       assert(len(data) == 2)
#       self.addToBuf(self.ax, data[0])
#       self.addToBuf(self.ay, data[1])
#
#   # update plot
#   def update(self, frameNum, a0, a1):
#       try:
#           line = self.ser.readline()
#           data = [float(val) for val in line.split()]
#           # print data
#           if(len(data) == 2):
#               self.add(data)
#               a0.set_data(range(self.maxLen), self.ax)
#               a1.set_data(range(self.maxLen), self.ay)
#       except KeyboardInterrupt:
#           print('exiting')
#
#       return a0,
#
#   # clean up
#   def close(self):
#       # close serial
#       self.ser.flush()
#       self.ser.close()
#
# # main() function
# def main():
#   # create parser
#   parser = argparse.ArgumentParser(description="LDR serial")
#   # add expected arguments
#   parser.add_argument('--port', dest='port', required=True)
#
#   # parse args
#   args = parser.parse_args()
#
#   #strPort = '/dev/tty.usbserial-A7006Yqh'
#   strPort = args.port
#
#   print('reading from serial port %s...' % strPort)
#
#   # plot parameters
#   analogPlot = AnalogPlot(strPort, 100)
#
#   print('plotting data...')
#
#   # set up animation
#   fig = plt.figure()
#   ax = plt.axes(xlim=(0, 100), ylim=(0, 1023))
#   a0, = ax.plot([], [])
#   a1, = ax.plot([], [])
#   anim = animation.FuncAnimation(fig, analogPlot.update,
#                                  fargs=(a0, a1),
#                                  interval=50)
#
#   # show plot
#   plt.show()
#
#   # clean up
#   analogPlot.close()
#
#   print('exiting.')
#
#
# # call main
# if __name__ == '__main__':
#     main()



import serial #Import Serial Library
import matplotlib.pyplot as plt
import numpy as np

arduinoSerialData = serial.Serial('/dev/ttyS8',115200) #Create Serial port object called arduinoSerialData

def ppause(event):
    pass

n=0
Datalist = np.zeros(5000)
fig = plt.figure()
ax = fig.add_subplot(111)
line, = ax.plot(Datalist)
ax.set_ylim(0,5000)

axpause = plt.axes([0.7, 0.05, 0.1, 0.075])
bpause = plt.Button(axpause, 'Pause')
bpause.on_clicked(ppause)

while (n<=10000):
    if (arduinoSerialData.inWaiting()>0):
        myData = arduinoSerialData.readline()
        # print(myData)
        try:
            A1 = myData.decode('utf-8').rstrip('\r\n').split(" ")[1]
            Datalist = np.append(Datalist[1:], float(A1))
            line.set_ydata(Datalist)
        except:
            pass
        plt.pause(0.05)
    n=n+1
