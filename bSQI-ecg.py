import matplotlib.pyplot as plt
import numpy as np
import pandas as pd 
import argparse
import neurokit2 as nk2
from datetime import datetime, timedelta

def parse_opts():
    parser = argparse.ArgumentParser() #Create parser
    parser.add_argument('-f', '--file',
                    type=str,
                    help="Path to the ECG file of a patient")
    parser.add_argument('-m1', '--method1',
                    default='hamilton',
                    type=str,
                    help="Method1 for R-peak detection")
    parser.add_argument('-m2', '--method2',
                    default='rodrigues',
                    type=str,
                    help="Method2 for R-peak detection") 
    args = parser.parse_args() # Parse arguments
    return args

class  ECG:
    def __init__(self, path):
        self.path = path
        self.rawdata = pd.read_csv(path)
        self.timestamps = self.rawdata.t
        self.voltage = self.rawdata.voltage

#HELPER -------------------------------------------------------------------------------

    #10s windows -------------------------------------------------------------------------
    # converts list of unix timestamp to list of datetime timestamps
    def timestampsCSV(self):
        times = self.voltage.tolist()
        res = [datetime.utcfromtimestamp(times[i]/1000) for i in range(len(times))]
        return res

    # 10s sized windows: returns the index of the timestamps (includes all indices from 0 to len(data))
    # list of relative timestamps -> list of indices of the windows [[10s window] [10s window] [10s window]...]
    def smallsamplesCSV(lst):
        start = 0
        end = 10
        temprange = np.arange(start,end,0.001)
        roundtemprange = np.round(temprange,3)  #nochmal np.round, da Probleme mit arange: hat Zahlen mit mehr Nachkommastellen generiert
        res = []
        liste = []
        for i in range(len(lst)):
            if lst[i] in roundtemprange:
                liste.append(i)
            else:
                res.append(liste)
                liste = []
                start += 10
                end += 10
                temprange = np.arange(start,end,0.001)
                roundtemprange = np.round(temprange,3)
                liste.append(i)
        res.append(liste)
        return res

    # 10s sized windows for CSV
    # data -> list of indices [[10s window] [10s window] [10s window]...]
    def datasamples(self):
        times = self.timestamps.to_numpy()
        timestamps = np.round(((times-times[0])/1000),3)
        res =  ECG.smallsamplesCSV(timestamps)
        if timestamps[max(res[-1])]-timestamps[min(res[-1])] < 8: # wenn das letzte sample kleiner als 8s ist, dann rausschmeißen und keinen bSQI berechnen
            res.pop()
        return res

    # Gaps --------------------------------------------------------------------------------
    # returns list of time differences in seconds between 2 r-peaks one method detected in ECG
    def difference_ts(method):
        rpeaks_ts = method[0]
        diff_lst = []
        i = 0
        while i < len(rpeaks_ts)-1:
            diff = (rpeaks_ts[i+1]-rpeaks_ts[i])
            diff_lst.append(diff)
            i += 1
        res = [diff_lst[j].total_seconds() for j in range(len(diff_lst))]
        return res

    # computes moving average for the time differences between r-peaks one method detects to get the gaps sizes between r-peaks
    def moving_avg(method):
        res, i, isnan = 0, 0, 1
        diff = ECG.difference_ts(method)
        df = pd.DataFrame(diff,columns=['sec'])
        res = ((df.rolling(11).mean()))["sec"].tolist() 
        while isnan:    # ersten 10 Werte setzen 
            if np.isnan(res[i]):
                res[i] = 1.5   
                i += 1
            else:
                isnan = 0
        return res

    # returns relative timestamps of the starting r-peak of the gap, if it is > moving average * 1.5, one method only
    # 1s = 1.000.000 microsec
    def rpeak_gap(self,method):
        diff = ECG.difference_ts(method)
        movingavg =  ECG.moving_avg(method)
        res = []
        for i in range(len(diff)):
            if diff[i] > (movingavg[i] * 1.5): 
                temp = method[1][i]
                res.append(np.round((self.timestamps[temp]-self.timestamps[0])/1000,2))
        return res  

    # returns relative timestamps of the gaps that both methods detect
    # rpeak_gap(), rpeak_gap() -> list of relative timestamps, where the gaps are in data (start of the gap)
    def gaps(gap1,gap2):
        res = []
        if len(gap1) > len(gap2):
            for i in range(len(gap1)):
                for j in range(len(gap2)):
                    if np.round(gap1[i],0) == np.round(gap2[j],0):
                        res.append(gap1[i])
                    elif gap1[i]-gap2[j] <= 0.15 and gap1[i]-gap2[j] > 0: 
                        res.append(gap1[i])
        else: 
            for i in range(len(gap2)):
                for j in range(len(gap1)):
                    if np.round(gap2[i],0) == np.round(gap1[j],0):
                        res.append(gap2[i])
                    elif gap2[i]-gap1[j] <= 0.15 and gap2[i]-gap1[j] > 0:
                        res.append(gap2[i])
        return res

    # returns list of indices of the 10s windows that have gaps
    # data, gaps() -> list 
    def samples_gap(self,gaps):
        start = 0
        stop = 11
        ranges = []
        res = []
        length = len(ECG.datasamples(self))
        while length != 0:
            ranges.append(np.round(np.arange(start,stop,0.01),2))
            start += 10
            stop += 10
            length -= 1
        for i in range(len(gaps)):
            for j in range(len(ranges)):
                if gaps[i] in ranges[j] and min(ranges[j]) not in res:
                    res.append(min(ranges[j])/10)
                    break
        return res

    # R-peaks -------------------------------------------------------------------------------
    # datetime timestamps of r-peaks detected by a method in every 10s window
    # method() + datasamples() -> list of datetime timestamps [[10s window] [10s window] [10s window]...]
    def rpeakstsinsample(method,samps): 
        rpeaks = method[1].copy() # Liste mit indices [...]
        timestamps = method[0].copy() # Liste mit timestamps [...]
        res = []
        empty = []
        for i in range(len(samps)):
            liste = []
            while len(rpeaks) != 0:
                if rpeaks[0] in samps[i]:
                    liste.append(timestamps[0])
                    rpeaks.pop(0)
                    timestamps.pop(0)
                else: 
                    res.append(liste)
                    break
            if len(rpeaks) == 0 and len(res) < len(samps): # keine r-peaks mehr gefunden aber die messung geht noch weiter -> leere Liste hinzufügen
                res.append(empty)
        return res

    # indices of r-peaks found by a method in every 10s window
    # method(), datasamples() -> list of indices of the r-peaks [[10s window] [10s window] [10s window]...]
    def rpeaksidxinsample(method,samps):  
        rpeaks = method[1].copy()
        res = []
        for i in range(len(samps)):
            liste = []
            while len(rpeaks) != 0:
                if rpeaks[0] in samps[i]:
                    liste.append(rpeaks[0])
                    rpeaks.pop(0)
                else: 
                    res.append(liste)
                    break
        return res

#bSQI -------------------------------------------------------------------------------

    # counter for number of r-peaks that both methods detect in each 10s-sample
    # data, method(), method(), datasamples() -> list of ints (res for N_matched)
    def N_matched(self, method1, method2, samps): 
        rpeaksmethod1 = ECG.rpeakstsinsample(method1, samps)
        rpeaksmethod2 = ECG.rpeakstsinsample(method2, samps)
        nmatched = []
        interval = timedelta(milliseconds=150) # max. 150ms-gap between two r-peaks
        for i in range(len(rpeaksmethod1)):
            ctr = 0
            if len(rpeaksmethod1[i]) > len(rpeaksmethod2[i]):
                used = []   #Prüfung, dass bspw. method1 für r-peak von method2 2 matches hat |I| (oder umgekehrt)
                for j in range(len(rpeaksmethod1[i])):
                    for l in range(len(rpeaksmethod2[i])):
                        if l in used:
                            continue
                        if rpeaksmethod1[i][j] > rpeaksmethod2[i][l] and ((rpeaksmethod1[i][j] - rpeaksmethod2[i][l]) <= interval) or (rpeaksmethod1[i][j] < rpeaksmethod2[i][l] and ((rpeaksmethod2[i][l] - rpeaksmethod1[i][j]) <= interval)) or (rpeaksmethod1[i][j] == rpeaksmethod2[i][l]): 
                            ctr += 1
                            used.append(l)
                            break
                nmatched.append(ctr)
            else: 
                used = []
                for j in range(len(rpeaksmethod2[i])):
                    for l in range(len(rpeaksmethod1[i])):
                        if l in used:
                            continue
                        if rpeaksmethod2[i][j] > rpeaksmethod1[i][l] and ((rpeaksmethod2[i][j] - rpeaksmethod1[i][l]) <= interval) or (rpeaksmethod2[i][j] < rpeaksmethod1[i][l] and ((rpeaksmethod1[i][l] - rpeaksmethod2[i][j]) <= interval)) or (rpeaksmethod1[i][l] == rpeaksmethod2[i][j]):
                            ctr += 1
                            used.append(l)
                            break
                nmatched.append(ctr)
        return nmatched

    # counter for all detected beats without double counting the matched ones: N_all = N_method1 + N_method2 - N_matched
    # data, method(), method(), datasamples() -> list of ints (res for N_all)
    def N_all(self,method1, method2,samps): 
        nmatched = ECG.N_matched(self,method1,method2,samps)
        rpeaksmethod1 = ECG.rpeakstsinsample(method1,samps)
        method1_rpeaks = [len(rpeaksmethod1[i]) for i in range(len(rpeaksmethod1))]
        rpeaksmethod2 = ECG.rpeakstsinsample(method2,samps)
        method2_rpeaks = [len(rpeaksmethod2[i]) for i in range(len(rpeaksmethod2))]
        res = [(method1_rpeaks[i] + method2_rpeaks[i] - nmatched[i]) for i in range(len(samps))]
        return res

    # computation bSQI = N_matched/N_all
    # data, method(), method() -> list of ints: bSQI values for every 10s window
    def bsqi(self,method1,method2): 
        samples = ECG.datasamples(self)
        nmatched = ECG.N_matched(self,method1,method2,samples)
        nall = ECG.N_all(self,method1,method2,samples)
        res = []
        for i in range(len(samples)):
            if nall[i] != 0:
                temp = nmatched[i]/nall[i]
                res.append(temp)
            else:
                res.append(0)
        return np.round(res,4) 

# R-PEAK DETECTION METHODS-------------------------------------------------------------------------------

    # datetime timestamps and indices of r-peaks detected by method in measurement (not splitted into 10s windows)
    # data -> list with datetime timestamps and indices of R-peaks [[timestamps][indices]]
    def hamilton(self): 
        voltages = self.voltage.tolist()  
        _, rpeak = nk2.ecg_peaks(nk2.ecg_clean(voltages,sampling_rate=250,method="hamilton2002"), sampling_rate=250, method ="hamilton2002")
        rpeaks = list(rpeak["ECG_R_Peaks"])
        timestamps = [datetime.utcfromtimestamp(self.timestamps[rpeaks[i]]/1000) for i in range(len(rpeaks))]  
        return timestamps,rpeaks 

    def nabian(self):  
        voltages = self.voltage.tolist()
        _, rpeak = nk2.ecg_peaks(voltages, sampling_rate=250, method ="nabian2018")
        rpeaks = list(rpeak["ECG_R_Peaks"])
        timestamps = [datetime.utcfromtimestamp(self.timestamps[rpeaks[i]]/1000) for i in range(len(rpeaks))]  
        return timestamps,rpeaks 

    def rodrigues(self): 
        voltages = self.voltage.tolist()
        _, rpeak = nk2.ecg_peaks(voltages, sampling_rate=250, method ="rodrigues2021")
        rpeaks = list(rpeak["ECG_R_Peaks"])
        timestamps = [datetime.utcfromtimestamp(self.timestamps[rpeaks[i]]/1000) for i in range(len(rpeaks))]  
        return timestamps,rpeaks

# PLOTS -------------------------------------------------------------------------------

    # plots ECG signal (blue) with 10s window marker, bSQI values, r-peaks detected by method1 (red), r-peaks detected by method2 (green)
    # data, method(), method(), bSQI() -> plot
    def plotting(self,method1,method2,bsqi,r_peaks=True):
        fig, ax = plt.subplots(1, 1)
        plt.rcParams['font.size'] = '20'
        samples = ECG.datasamples(self)
        samples_ends = [max(samples[i]) for i in range(len(samples))]
        time = self.timestamps.to_numpy()
        timestamps = (time-time[0])/1000    #relative Zeit für x-Achse
        ax.plot(timestamps, self.voltage.to_numpy(), zorder=3)    # ECG
        ax.set_xlabel('Seconds', fontsize=25)
        ax.set_ylabel('Voltage', fontsize=25)
        if r_peaks:
            method1_timestamps = timestamps[method1[1]]
            method2_timestamps = timestamps[method2[1]]
            ends = timestamps[samples_ends]
            min_volt = self.rawdata.voltage.min()
            max_volt = self.rawdata.voltage.max()
            ax.vlines(ends, min_volt, max_volt,label='10 s-window marker',color='black',linestyles='dotted', zorder=1) #marker 10s
            ax.vlines(method1_timestamps, min_volt, max_volt,label='Method 1', color='red', zorder=1,linewidth=0.8) #r-peaks method1
            ax.vlines(method2_timestamps, min_volt, max_volt, label='Method 2', color='green', zorder=1,linewidth=0.8) #r-peaks method2
            ax.legend(loc='upper left') 
            for i, txt in enumerate(bsqi): 
                ax.annotate(txt, (ends[i], 0)) # bsqi values 
        plt.draw()           

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    args = parse_opts()
    ecg = ECG(args.file)
    method1, method2 = 0,0
    # R-peak detection for method1
    if args.method1 ==  'hamilton':
        method1 = ecg.hamilton() 
    elif args.method1 == 'rodrigues':
        method1 = ecg.rodrigues()
    else:
        method1 = ecg.nabian()
    # R-peak detection for method2
    if args.method2 ==  'hamilton':
        method2 = ecg.hamilton() 
    elif args.method2 == 'rodrigues':
        method2 = ecg.rodrigues()
    else:
        method2 = ecg.nabian()
    bsqi = ecg.bsqi(method1,method2)
    ecg.plotting(method1,method2,bsqi)
    plt.show()

