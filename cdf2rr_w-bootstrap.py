#!/usr/bin/python -u

'''
Recurrence Rate is modeled as the derivative of the mean cumulative number of events through time.
The mean cumulative number of events through time is found from several sets (Monte Carlo) of 
potential cumulative number of events through time.
'''

import numpy as np
import matplotlib.pyplot as plt
import sys

###################
#Use all MC solutions to model Cumulative Number of Events
def cumulative_function(data,time):
	'''
	Creates a mean value of cumulative events for all times in the time list.
	Each row in the data array should be a list of event times. Each row is
	treated as one Monte Carlo solution of event times, for a Volcanic Event
	Age Model (VEAM)
	'''
	rows = len(data)                #returns (length,width) of data
	simulation_bins = range(rows+1) #histogram bins for the cumulative functions
	CDF_func = np.zeros(len(time))  #Mean Cumulative Function list
	
	for i,t in enumerate(time):
		#Time Histogram...
		#    count number of events before time t in each data row
		histogram_t = np.histogram(np.where(data>t)[0],bins=simulation_bins)[0]
		#[row 1: 0 events, row 2: 5 events, row 3: 2 events, ...]
	
		#Average of events in MC solutions with time (CDF)
		
		CDF_func[i] = np.mean(histogram_t)
	
	return CDF_func


###################
#Use all MC solutions to model Cumulative Number of Volumes
def cumulative_volume(data,time,volumes):
	'''
	Creates a mean value of cumulative events for all times in the time list.
	Each row in the data array should be a list of event times. Each row is
	treated as one Monte Carlo solution of event times, for a Volcanic Event
	Age Model (VEAM)
	'''
	rows = len(data)                #returns (length,width) of data
	simulation_bins = range(rows+1) #histogram bins for the cumulative functions
	vol_cdf_func = np.zeros(len(time))  #Mean Cumulative Function list
	
	if len(np.shape(data))>1:
		for i,t in enumerate(time):
			#Time Histogram...
			#    count number of events before time t in each data row
			#histogram_t = np.histogram(np.where(data>t)[0],bins=simulation_bins)[0]
			#This might work:
			# np.sum(volumes[np.where(data>t)[1]]) gives total volume output of all MC sets combined.
			#np.sum(volumes[np.where(data>t)[1]])/rows gives mean volume output per MC set
		
			vol_cdf_func[i] = np.sum(volumes[np.where(data>t)[1]])/rows #gives mean volume output per MC set

			#[row 1: 0 events, row 2: 5 events, row 3: 2 events, ...]
	
			#Average of events in MC solutions with time (CDF)
		
			#CDF_func[t] = np.mean(histogram_t)
	
	else:
		for i,t in enumerate(time):
			vol_cdf_func[i] = np.sum(volumes[np.where(data>t)[0]])/rows #gives mean volume output per MC set
			
	return vol_cdf_func

###################
#Model RR by calculating derivative of mean Cumulative Distribution Function
def derivative(function):
	'''
	Calculates the derivative of a function using the Central difference method.
	Backward difference is used for the final element of the function,
	Forward difference is used for the first element.
	'''
	count = len(function)
	dfunction = np.zeros(count)
	for t in range(count):
		if t>0:
			if t<(count-1):
				dfunction[t] = (function[t-1]-function[t+1])/2 #Central  difference
			else:
				dfunction[t] = function[t-1]-function[t]       #Backward difference
		else:
			dfunction[t] = function[t]-function[t+1]         #Forward  difference
	return dfunction





#Get a list of events x 10k sets of potential ages
MCdata = np.loadtxt("crater_neighbor_1k.txt",skiprows=1,delimiter=",")
(MCsets,eventCount) = np.shape(MCdata)
print "Loaded MC dataset"

time = np.arange(350)           #year range (0-350 Ma, sample at each Ma)
tstp = 10
time5 = np.arange(0,350,tstp)
'''
Cumulative_events = cumulative_function(MCdata,time)
print "Found mean cumulative events function of all MC solutions"
Recurrence_rate   = derivative(Cumulative_events)
print "Found derivative of cumulative events with time"
'''

Cumulative_events5 = cumulative_function(MCdata,time5)
print "Found mean cumulative events function (5) of all MC solutions"
Recurrence_rate5   = derivative(Cumulative_events5)/tstp
print "Found derivative of cumulative events (5) with time"

#################RECURRENCE RATE UNCERTAINTY#######################
print "\n\n         ~~~~~~~~Event RR~~~~~~~~~"
CDF_each = np.zeros(( MCsets,len(time5) ))
for i in range(len(MCdata)):
	CDF_each[i] = cumulative_function(MCdata[i:(i+1)],time5)

#np.save("CDF_each",CDF_each)

#CDF_each = np.load("CDF_each5.npy")
print "Found cumulative event functions for each MC solution"

#"sawtooth" of each MC sol'n
RR_each = np.zeros(( MCsets,len(time5) ))
for i,s in enumerate(CDF_each):
	RR_each[i] = derivative(s)/tstp


#Compare 1Ma and 5Ma RR results
for i in range(1000):
	plt.plot(time5,RR_each[i],c='0.95')
plt.plot(time5,np.percentile(RR_each, 90, axis=0),c='r',label="95%")
plt.plot(time5,np.percentile(RR_each, 10, axis=0),c='r',label="5%")
plt.plot(time5,np.percentile(RR_each, 50, axis=0),c='g',label="50%")
#plt.plot(time,Recurrence_rate,c='k',label='CDF derivative RR Mean')
plt.plot(time5,Recurrence_rate5,c='k',label='CDF derivative RR Mean 5Ma')

#Chart attributes
plt.xlim([min(time),max(time)])
plt.gca().invert_xaxis()
plt.ylabel('Events per Myr')
plt.xlabel('Millions of years before present')
plt.title('Recurrence Rate Models of volcanism on Arsia Mons')
plt.legend(loc='upper left')
#plt.show()
plt.cla()

#################CUMULATIVE   UNCERTAINTY#######################
print "\n\n         ~~~~~~~~CDF of event count~~~~~~~~~"
#mean CDF of MC solutions
CDF_all = np.load("CDF_each.npy")

plt.plot(np.mean(CDF_all, axis=0),c='b',label='bootstrap CDF Mean')
plt.plot(np.percentile(CDF_all, 10, axis=0),c='b',ls='--',label='bootstrap CDF 10th per')
plt.plot(np.percentile(CDF_all, 90, axis=0),c='b',ls='--',label='bootstrap CDF 90th per')


#Chart attributes
plt.xlim([min(time),max(time)])
plt.gca().invert_xaxis()
plt.ylabel('Cumulative Events')
plt.xlabel('Millions of years before present')
plt.title('Cumulative Volcanism on Arsia Mons')
plt.legend(loc='upper left')
#plt.show()
plt.cla()

#################CUMULATIVE VOLUME UNCERTAINTY#####################
print "\n\n         ~~~~~~~~Incorporating Volumes~~~~~~~~~"
volumes = np.loadtxt("verrm_sorted_areas.txt",skiprows=1)
print volumes
print np.sum(volumes)
mean_cum_vol = cumulative_volume(MCdata,time5,volumes)

cum_vols = np.zeros(( MCsets,len(time5) ))
for i in range(MCsets):
	cum_vols[i] = cumulative_volume(MCdata[i],time5,volumes)
cum_vols*=len(volumes)

plt.plot(time5,mean_cum_vol,c='r',label='Vol CDF mean, at once')
plt.plot(time5,np.mean(cum_vols, axis=0),c='b',label='Vol CDF Mean')
plt.plot(time5,np.percentile(cum_vols, 10, axis=0),c='b',ls='--',label='Vol CDF 10th per')
plt.plot(time5,np.percentile(cum_vols, 90, axis=0),c='b',ls='--',label='Vol CDF 90th per')

#Chart attributes
plt.xlim([min(time),max(time)])
plt.gca().invert_xaxis()
plt.ylabel('Cumulative Volume')
plt.xlabel('Millions of years before present')
plt.title('Cumulative Volcanism on Arsia Mons')
plt.legend(loc='upper left')
plt.show()
plt.cla()

#################VOLUME FLUX UNCERTAINTY#####################
print "\n\n         ~~~~~~~~Volume Flux over time!~~~~~~~~~"

#cum_vols
#np.save("CDF_each",CDF_each)

#CDF_each = np.load("CDF_each5.npy")
#print "Found cumulative event functions for each MC solution"

#"sawtooth" of each MC sol'n
VF_each = np.zeros(( MCsets,len(time5) ))
for i,s in enumerate(cum_vols):
	VF_each[i] = derivative(s)/tstp


#Compare 1Ma and 5Ma RR results
for i in range(1000):
	plt.plot(time5,VF_each[i],c='0.95')
plt.plot(time5,np.percentile(VF_each, 90, axis=0),c='r',label="95%")
plt.plot(time5,np.percentile(VF_each, 10, axis=0),c='r',label="5%")
plt.plot(time5,np.percentile(VF_each, 50, axis=0),c='g',label="50%")
#plt.plot(time,Recurrence_rate,c='k',label='CDF derivative RR Mean')
plt.plot(time5,np.mean(VF_each, axis=0),c='k',label='CDF derivative VF Mean 5Ma')

#Chart attributes
plt.xlim([min(time),max(time)])
plt.gca().invert_xaxis()
plt.ylabel('Events per Myr')
plt.xlabel('Millions of years before present')
plt.title('Recurrence Rate Models of volcanism on Arsia Mons')
plt.legend(loc='upper left')
plt.show()
plt.cla()










