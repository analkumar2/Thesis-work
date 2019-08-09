V_trace = []
time = []
with open('300pA_B2012') as file:
    reader = csv.reader(file, delimiter = '\t')
    a = next(reader)
    a = next(reader)
    for row in reader:
        time.append(float(row[0]))
        V_trace.append(float(row[1]))
        
plt.plot(time, V_trace)
plt.xlabel('Time (ms)')
plt.ylabel('Membrane potential (mV)')
plt.title('300pA current injection Bianchi2012 full model')
plt.show()
