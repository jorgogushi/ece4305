import time
import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy import signal
import sys as sys
import scipy.integrate as integrate
from scipy.signal import butter,lfilter, freqz, filtfilt

data_array = np.loadtxt("/Users/jorgogushi/Desktop/data_array_1fs_2426fc_1", dtype=np.cdouble)


sample_rate = 1e6 #Hz
fc = 2426e6 #Hz
buffer_size = len(data_array)
#start = int(8e4)
#end = int(8e4+buffer_size)
#print(start)
#print(end)
#data_array = data_array[start:end]
#data_array = buffer_size

f, t, Sxx = signal.spectrogram(data_array, sample_rate, return_onesided=False)
f = np.fft.fftshift(f)+fc
Sxx = np.fft.fftshift(Sxx, axes=0,)
Sxx = np.transpose(Sxx)
Sxx = np.flipud(Sxx)

out_data = [np.abs(x)*np.sign(np.angle(x)) for x in data_array]
time_domain = np.linspace(0,buffer_size/sample_rate,len(data_array))
freq_domain = np.linspace(fc-sample_rate/2,fc+sample_rate/2,buffer_size)

shifted_fft= np.abs(np.fft.fftshift(np.fft.fft((data_array))))

plt.subplot(2,1,1)
plt.plot(time_domain, out_data)
plt.xlabel("Time [sec]")
plt.ylabel("Magnitude")

plt.subplot(2,1,2)
plt.pcolormesh(f, t, Sxx, shading="gouraud")
plt.xlabel("Freqeuency [Hz]")
plt.ylabel("Time [sec]")
plt.show()


#Coarse Frequency Correction
f_c = int(0.5*len(freq_domain))

half_int_stop_point = 0.5*np.sum(shifted_fft)
current_rolling_integral = 0
for i in range(len(shifted_fft)):
    current_rolling_integral += shifted_fft[i]
    if current_rolling_integral > half_int_stop_point:
        #calculate f_offset based on i here
        freq_offset = freq_domain[i] - fc
        break

samples_shifted = data_array*np.exp(1j*2*np.pi*freq_offset*time_domain)

print(freq_offset)
samples_of_f_1 = np.abs(np.fft.fftshift(np.fft.fft((samples_shifted))))
fig, (plotT, plotF) = plt.subplots(2)
plotT.plot(freq_domain,  samples_of_f_1)
plotF.plot(freq_domain, shifted_fft)
plt.show()

def RunningAVG(data,offset):
    sum =0
    for i in range(100):
        sum += np.abs(data[i+offset])
    return sum/100

startPoint= 0
endPoint = 0
difference = 0
ran=False
for i in range(len(samples_shifted)-100):
    #print(i)
    #print(RunningAVG(out_data,i))
    if (RunningAVG(samples_shifted,100)>600 and ran==False):
        print("entered1")
        j=1
        while(RunningAVG(samples_shifted,100+j)>600):
            j+=1
        i=j
        ran = True
    else:
        print("entered")
        if (RunningAVG(samples_shifted,i)>600):
            startPoint = i-100
            g = 1
            while(RunningAVG(samples_shifted,i+g)>600):
                g+=1
            endpoint = i+g+100
            difference = g+300
            print("entered1")
            break
packetdata = []
print(difference)
print(startPoint)
print(endPoint)
#print(len(data_array))
for i in range(difference):
    print(i)
    packetdata.append(samples_shifted[startPoint+i])
    

out_packetdata = [np.abs(x)*np.sign(np.angle(x)) for x in packetdata]
new_time_domain = np.linspace(0,difference/sample_rate,len(packetdata))
plt.plot(new_time_domain, out_packetdata)
plt.xlabel("Time [sec]")
plt.ylabel("Magnitude")
plt.show()

# DIGITAL PHASE LOCK LOOP #

# Phase Error Detector (PED)

samples_of_dpll = np.zeros(len(packetdata)) 
phase = np.angle(packetdata) 

samples_of_dpll_angle = np.angle(packetdata)
samples_of_dpll_magnitude = np.abs(packetdata)
real_part = np.zeros(len(packetdata))
imaginary_part = np.zeros(len(packetdata))

for j in range(len(packetdata)):
    if -np.pi < phase[j] < -np.pi/2:
        phase[j] = -np.pi - phase[j]
    elif np.pi/2 < phase[j] < np.pi:
        phase[j] = np.pi - phase[j]

    avg = phase[j-1]

    if -np.pi < samples_of_dpll_angle[j] < -np.pi/2:
        samples_of_dpll_angle[j] = samples_of_dpll_angle[j] + avg
    elif np.pi/2 < samples_of_dpll_angle[j] < np.pi:
        samples_of_dpll_angle[j] = samples_of_dpll_angle[j] + avg
    else:
        samples_of_dpll_angle[j] = samples_of_dpll_angle[j] - avg

    real_part[j] = samples_of_dpll_magnitude[j]*np.cos(samples_of_dpll_angle[j])
    imaginary_part[j] = samples_of_dpll_magnitude[j]*np.sin(samples_of_dpll_angle[j])

samples_of_dpll = real_part + 1j*imaginary_part

plt.scatter(np.real(packetdata),np.imag(packetdata),marker = '.')
plt.scatter(real_part, imaginary_part, marker = 'x', color = 'r')
plt.show()

# Converting to Binary Data
positive_phase = np.zeros(len(samples_of_dpll))

for i in range(len(samples_of_dpll)):
    for j in range(len(samples_of_dpll)):
        positive_phase[j]=np.angle(packetdata[j])
        if positive_phase[i]<0:
            positive_phase[j] = 2*np.pi + positive_phase[j]

# Frequency Shifting

angle_change = np.zeros(len(samples_of_dpll))

angle_change[0] = positive_phase[0]
for i in range(1, len(samples_of_dpll)):
    angle_change[i] = positive_phase[i] - positive_phase[i-1]


# Actual Conversion

binary_samples = np.zeros(len(angle_change))

for i in range(len(binary_samples)):
    if angle_change[i] > 0:
        binary_samples[i] = 1
    if angle_change[i] < 0:
        binary_samples[i] = 0

print(binary_samples)

# Preamble Logic Implementation

array = np.array(binary_samples)
for k in range(len(array)-7) :
    if array[k] == 0 and \
    array[k+1] == 1 and \
    array[k+2] == 0 and \
    array[k+3] == 1 and \
    array[k+4] == 0 and \
    array[k+5] == 1 and \
    array[k+6] == 0 and \
    array[k+7] == 1 or \
    array[k] == 1 and \
    array[k+1] == 0 and \
    array[k+2] == 1 and \
    array[k+3] == 0 and \
    array[k+4] == 1 and \
    array[k+5] == 0 and \
    array[k+6] == 1 and \
    array[k+7] == 0: 
        find_preamble = k
        break
print('The index where the Preamble starts is:',find_preamble)

access_adress = array[48:80]
print('The Access Adress is:',access_adress)

# Dewhitening

def str_xor(a,b): #returns a list of bits
	return list(map(lambda x: 0 if x[0] is x[1] else 1, zip(a,b)))

def bit_xor(a,b): #returns a list of bits
	return list(map(lambda x: x[0] ^ x[1], zip(a,b)))

def dewhiten_bits(bits, channel_num):
	front_half = [1,1,0,0]
	back_half = [1,1,0]
	if channel_num == 37:
		back_half = [1,0,1]
	elif channel_num == 38:
		back_half = [1,1,0]
	elif channel_num == 39:
		back_half = [1,1,1]
	else:
		print("you didn't call this correctly")
		quit()
	#LSB on left, initialize to [1, channel in binary]
	current_state = [front_half,back_half] #output of lfsr on right
	lfsr_out_str = ""
	lfsr_out_bit = []
	for i in range(len(bits)):
		out_bit = current_state[1][-1]
		lfsr_out_str = lfsr_out_str + str(out_bit)
		lfsr_out_bit.append(out_bit)
		current_state[1] = [current_state[0][-1] ^ out_bit] + current_state[1][:-1]
		current_state[0] = [out_bit] + current_state[0][:-1]
	return str_xor(bits, lfsr_out_str)
	return bit_xor(bits, lfsr_out_bit)