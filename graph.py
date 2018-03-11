import matplotlib.pyplot as plt
import re

data_file = open("file_oper_refuse.txt")
data_file2 = open("file_load_koef.txt")
data_file3 = open("file_reverse.txt")
data_file4 = open("file_load_koef2.txt")
data_file5 = open("file_len_queue.txt")
data_file6 = open("file_time_queue.txt")
data_file7 = open("file_load_koef5.txt")
data_file8 = open("file_len_queue5.txt")
data_file9 = open("file_time_queue5.txt")

X = []
Y = []
X1 = []
Y1 = []
X2 =[]
Y2 = []
X3 = []
Y3 = []
X4 = []
Y4 = []
X5 = []
Y5 = []

X6 = []
Y6 = []
X7 = []
Y7 = []
X8 = []
Y8 = []

for r in data_file:
        res = re.findall(r"[-+]?\d*\.\d+|\d+", r)
        #print res
        X.append(float(res[0]))
        Y.append(float(res[1]))

for r in data_file2:
        res = re.findall(r"[-+]?\d*\.\d+|\d+", r)
        #print res
        X1.append(float(res[0]))
        Y1.append(float(res[1]))

for r in data_file3:
        res = re.findall(r"[-+]?\d*\.\d+|\d+", r)
        print res
        X2.append(float(res[0]))
        Y2.append(float(res[1]))

for r in data_file4:
        res = re.findall(r"[-+]?\d*\.\d+|\d+", r)
        print res
        X3.append(float(res[0]))
        Y3.append(float(res[1]))

for r in data_file5:
        res = re.findall(r"[-+]?\d*\.\d+|\d+", r)
        print res
        X4.append(float(res[0]))
        Y4.append(float(res[1]))

for r in data_file6:
        res = re.findall(r"[-+]?\d*\.\d+|\d+", r)
        print res
        X5.append(float(res[0]))
        Y5.append(float(res[1]))

for r in data_file7:
        res = re.findall(r"[-+]?\d*\.\d+|\d+", r)
        print res
        X6.append(float(res[0]))
        Y6.append(float(res[1]))

for r in data_file8:
        res = re.findall(r"[-+]?\d*\.\d+|\d+", r)
        print res
        X7.append(float(res[0]))
        Y7.append(float(res[1]))

for r in data_file9:
        res = re.findall(r"[-+]?\d*\.\d+|\d+", r)
        print res
        X8.append(float(res[0]))
        Y8.append(float(res[1]))


plt.plot(X, Y)
plt.show()
plt.plot(X1, Y1)
plt.show()
plt.plot(X2, Y2)
plt.show()
plt.plot(X3, Y3)
plt.show()
plt.plot(X4, Y4)
plt.show()
plt.plot(X5, Y5)
plt.show()
plt.plot(X6, Y6)
plt.show()
plt.plot(X7, Y7)
plt.show()
plt.plot(X8, Y8)
plt.show()
