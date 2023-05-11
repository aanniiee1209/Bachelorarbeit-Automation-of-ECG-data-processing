import matplotlib.pyplot as plt
import numpy as np
import pandas as pd 
import os
import fnmatch
import argparse
import neurokit2 as nk2
from datetime import datetime, timedelta
import seaborn as sns

def parse_opts():
    parser = argparse.ArgumentParser() #Create parser
    parser.add_argument('-f', '--file',
                    type=str,
                    help="Path to the ECG files of a patient")
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
    # scatter plot for bSQI values (each combination of methods) and markers for gaps in the data
    # bsqi(), samples_gap(), bsqi(), samples_gap(), bsqi(), samples_gap() -> plot
    def plot_bsqis(bsqi1,gaps1,bsqi2,gaps2,bsqi3,gaps3):
        fig, (ax1, ax2,ax3) = plt.subplots(3, 1,figsize=(20, 11))
        plt.rcParams['font.size'] = '40'
        gaps1 = []
        for i in range(len(gaps1)):
            gaps1.append(int(gaps1[i]))
        for j in range(len(gaps1)):
            ax1.scatter(gaps1[j],bsqi1[gaps1[j]], marker = 's', zorder = 10, facecolors='none', edgecolors ='black', s=200)
        ax1.plot(bsqi1,'o--', color='red', label='Hamilton & Rodrigues')
        ax1.set_xticks(np.arange(0, len(bsqi1),10))
        ax1.set_yticks(np.arange(0, 1.1 ,0.25))   
        ax1.tick_params(axis='both', labelsize=20)
        ax1.grid()
        ax1.set_title('Hamilton & Rodrigues', fontsize=30)
        gaps2 = []
        for i in range(len(gaps2)):
            gaps2.append(int(gaps2[i]))
        for j in range(len(gaps2)):
            ax2.scatter(gaps2[j],bsqi2[gaps2[j]], marker = 's', zorder = 10, facecolors='none', edgecolors ='black', s=200)
        ax2.plot(bsqi2,'o--', color='green', label='Hamilton & Nabian')
        ax2.set_ylabel('bSQI', fontsize=35)
        ax2.set_xticks(np.arange(0, len(bsqi1),10))
        ax2.set_yticks(np.arange(0, 1.1 ,0.25)) 
        ax2.tick_params(axis='both', labelsize=20) 
        ax2.grid() 
        ax2.set_title('Hamilton & Nabian', fontsize=30)
        gaps3 = []
        for i in range(len(gaps3)):
            gaps3.append(int(gaps3[i]))
        for j in range(len(gaps3)):
            ax3.scatter(gaps3[j],bsqi3[gaps3[j]], marker = 's', zorder = 10, facecolors='none', edgecolors ='black', s=200)
        ax3.plot(bsqi3,'o--', color='blue', label='Rodrigues & Nabian')
        ax3.set_xlabel('10 s-windows', fontsize=35)
        ax3.set_xticks(np.arange(0, len(bsqi1),10))
        ax3.set_yticks(np.arange(0, 1.1 ,0.25))
        ax3.tick_params(axis='both', labelsize=20)
        ax3.grid()   
        ax3.set_title('Rodrigues & Nabian', fontsize=30)
        plt.tight_layout(pad = 1.8)
    
    # plots histogram for bSQI values (each combination of methods)
    # bsqi(), bsqi(), bsqi() -> plot
    def plot_histos(bsqi1,bsqi2,bsqi3):
        fig, ax = plt.subplots(1, 1,figsize=(20, 11))
        plt.rcParams['font.size'] = '30'
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontsize(25)
        temp = [] 
        temp.append(bsqi1)    
        temp.append(bsqi2) 
        temp.append(bsqi3) 
        colors = ['red','green','blue']
        label= ['Hamilton & Rodrigues','Hamilton & Nabian','Rodrigues & Nabian']
        ax.hist(temp, bins= [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1], color=colors,label=label, ec='black')
        ax.set_axisbelow(True)
        ax.grid(axis='y')
        plt.xlabel('bSQI', fontsize=30)
        plt.ylabel('#10 s-windows',fontsize=30)
        plt.xticks(np.arange(0,1.1,0.1))
        plt.yticks(np.arange(0,len(bsqi1),25))
        plt.legend() 
        plt.tight_layout() 

    # plot heatmap for sorted bSQI values (each combination of methods)
    # bsqi(), bsqi(), bsqi() -> plot
    def plot_heatmapssorted(bsqi1,bsqi2,bsqi3):
        bsqis = [] 
        bsqis.append(sorted(bsqi1))    
        bsqis.append(sorted(bsqi2))
        bsqis.append(sorted(bsqi3))
        fig, ax = plt.subplots(figsize=(14, 7))
        plt.rcParams['font.size'] = '20'
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontsize(16)
        ylabels = ['a)', 'b)', 'c)'] #a) hamilton_rodrigues b) hamilton_nabian c) rodrigues_nabian 
        hm = sns.heatmap(bsqis,cbar_kws={'ticks': [0.0, 0.25, 0.5, 0.75, 1.0]}, vmin=0, vmax=1, xticklabels = 25, yticklabels=ylabels)
        hm.set_yticklabels(labels=hm.get_yticklabels(), va='center')
        plt.title('Heatmap - Sorted bSQIs')
        plt.yticks(rotation=0) 

    # plot heatmap for bSQI values (not sorted, each combination of methods)
    # bsqi(), bsqi(), bsqi() -> plot
    def plot_heatmaps(bsqi1,bsqi2,bsqi3):
        bsqis = [] 
        bsqis.append(bsqi1)    
        bsqis.append(bsqi2) 
        bsqis.append(bsqi3)
        fig, ax = plt.subplots(figsize=(14, 7))
        plt.rcParams['font.size'] = '20'
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontsize(16)
        ylabels = ['a)', 'b)', 'c)'] #a) hamilton_rodrigues b) hamilton_nabian c) rodrigues_nabian    
        hm = sns.heatmap(bsqis,cbar_kws={'ticks': [0.0, 0.25, 0.5, 0.75, 1.0]}, vmin=0, vmax=1, xticklabels = 25, yticklabels=ylabels)
        hm.set_yticklabels(labels = hm.get_yticklabels(), va='center')
        plt.title('Heatmap - bSQIs')
        plt.yticks(rotation=0) 

    # plot heatmap for amount of bSQI values, that have the same range
    # bsqi(), bsqi(), bsqi() -> plot
    def plot_bsqisheatmap(bsqi1,bsqi2,bsqi3): 
        bsqis = []
        val1 = plt.hist(bsqi1,bins = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])[0]
        bsqis.append(val1)
        val2 = plt.hist(bsqi2,bins = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])[0]
        bsqis.append(val2)
        val3 = plt.hist(bsqi3,bins = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])[0]
        bsqis.append(val3)
        fig, ax = plt.subplots(figsize=(14, 7))
        plt.rcParams['font.size'] = '20'
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontsize(16)
        xlabels = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
        ylabels = ['a)', 'b)', 'c)'] #a) hamilton_rodrigues b) hamilton_nabian c) rodrigues_nabian  
        hm = sns.heatmap(bsqis,xticklabels = xlabels, yticklabels=ylabels,annot=True, fmt='.1f')
        hm.set_yticklabels(labels=hm.get_yticklabels(), va='center')
        plt.title('Heatmap')
        plt.yticks(rotation=0) 

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    args = parse_opts()
    flist = fnmatch.filter(os.listdir(str(args.file)), "*.csv") # data files in current directory
    if not(os.path.exists(args.file + '/bsqi')):    # create new directory for results 
        os.mkdir(args.file + '/bsqi')
    currentdir = os.getcwd()    # current working directory
    datafiles = len(flist)
# STATS variables for all data files of one patient -----------------------------------------------
    samples_total = 0 
    length_total = 0
    # R-peak detection methods
    hamilton_total = 0
    rodrigues_total = 0
    nabian_total = 0
    #Hamilton & Rodrigues             
    hamiltonrodrigues_avgtotal = []  
    hrcorrupted_total = 0
    hrgut_total = 0
    hrmittel_total = 0
    hrschlecht_total = 0
    #Hamilton & Nabian             
    hamiltonnabian_avgtotal = []
    hncorrupted_total = 0
    hngut_total = 0
    hnmittel_total = 0
    hnschlecht_total = 0
    #Rordrigues & Nabian            
    rodriguesnabian_avgtotal = []
    rncorrupted_total = 0
    rngut_total = 0
    rnmittel_total = 0
    rnschlecht_total = 0
    for i in range(datafiles):
# STATS for each data file of a patient
        data = args.file + "/" + str(flist[i])  # file path 
        ecg =  ECG(data)    # ECG
        length = (ecg.timestamps[len(ecg.rawdata)-1] - ecg.timestamps[0])/1000  # data length
        length_total += length
        samples = ecg.datasamples() 
        samples_total += len(samples) 
        #Hamilton & Rodrigues stats            
        hrgut = 0
        hrmittel = 0
        hrschlecht = 0
        #Hamilton & Nabian stats        
        hngut = 0
        hnmittel = 0
        hnschlecht = 0
        #Rordrigues & Nabian stats          
        rngut = 0
        rnmittel = 0
        rnschlecht = 0
        # R-PEAK DETECTION -------------------------------------------------
        hamilton = ecg.hamilton() 
        hamilton_total += len(hamilton[1]) 
        rodrigues = ecg.rodrigues()
        rodrigues_total += len(rodrigues[1])
        nabian = ecg.nabian()
        nabian_total += len(nabian[1])
        # bSQI -------------------------------------------------------------
        hamilton_rodrigues = ecg.bsqi(hamilton,rodrigues)
        hamilton_nabian = ecg.bsqi(hamilton,nabian)
        rodrigues_nabian = ecg.bsqi(rodrigues,nabian)
        bsqis = []  # all bSQI values 
        bsqis.append(hamilton_rodrigues)
        bsqis.append(hamilton_nabian)
        bsqis.append(rodrigues_nabian)
        # GAPS -------------------------------------------------------------
        hamilton_gaps = ecg.rpeak_gap(hamilton)
        rodrigues_gaps = ecg.rpeak_gap(rodrigues)
        nabian_gaps = ecg.rpeak_gap(nabian)
        hrgaps =  ECG.gaps(hamilton_gaps,rodrigues_gaps)
        hngaps =  ECG.gaps(hamilton_gaps,nabian_gaps)
        rngaps =  ECG.gaps(rodrigues_gaps,nabian_gaps)
        hrgapsamp = ecg.samples_gap(hrgaps)
        hngapsamp = ecg.samples_gap(hngaps)
        rngapsamp = ecg.samples_gap(rngaps)
        #percentage of gaps in whole data
        hrgapscountpct = np.round((100/len(samples))*len(hrgaps),3) 
        hngapscountpct = np.round((100/len(samples))*len(hngaps),3)
        rngapscountpct = np.round((100/len(samples))*len(rngaps),3) 
        #total number of (possibly) corrupted samples (all data files)
        hrcorrupted_total += len(hrgaps)
        hncorrupted_total += len(hngaps)
        rncorrupted_total += len(rngaps)
        # bSQI STATS -----------------------------------------------
        #Hamilton & Rodrigues
        for j in range(len(hamilton_rodrigues)):
            if hamilton_rodrigues[j] >= 0.8:
                hrgut += 1
                hrgut_total += 1
            elif  0.8 > hamilton_rodrigues[j] >= 0.5:
                hrmittel += 1
                hrmittel_total += 1
            else:
                hrschlecht += 1
                hrschlecht_total += 1
        #Hamilton & Nabian
        for k in range(len(hamilton_nabian)):
            if hamilton_nabian[k] >= 0.8:
                hngut += 1
                hngut_total += 1
            elif  0.8 > hamilton_nabian[k] >= 0.5:
                hnmittel += 1
                hnmittel_total += 1
            else:
                hnschlecht += 1
                hnschlecht_total += 1
        #Rodrigues & Nabian
        for l in range(len(rodrigues_nabian)):
            if rodrigues_nabian[l] >= 0.8:
                rngut += 1
                rngut_total += 1
            elif  0.8 > rodrigues_nabian[l] >= 0.5:
                rnmittel += 1
                rnmittel_total += 1
            else: 
                rnschlecht += 1
                rnschlecht_total += 1
        hamiltonrodrigues_avgtotal.append(np.round(np.average(hamilton_rodrigues),3))
        hamiltonnabian_avgtotal.append(np.round(np.average(hamilton_nabian),3))
        rodriguesnabian_avgtotal.append(np.round(np.average(rodrigues_nabian),3))
# PLOTS -----------------------------------------------------
        ECG.plot_bsqis(hamilton_rodrigues, hrgapsamp, hamilton_nabian, hngapsamp, rodrigues_nabian, rngapsamp)
        plt.savefig(args.file + '/bsqi/' + str(flist[i]) + '_bsqis.png')
        plt.clf() 
        ECG.plot_histos(hamilton_rodrigues, hamilton_nabian, rodrigues_nabian)
        plt.savefig(args.file + '/bsqi/' + str(flist[i]) + '_histo.png')
        plt.clf() 
        ECG.plot_heatmapssorted(hamilton_rodrigues, hamilton_nabian, rodrigues_nabian)
        plt.savefig(args.file + '/bsqi/' + str(flist[i]) + '_heatmap-sorted.png') 
        plt.clf()
        ECG.plot_heatmaps(hamilton_rodrigues, hamilton_nabian, rodrigues_nabian)
        plt.savefig(args.file + '/bsqi/' + str(flist[i]) + '_heatmap.png') 
        plt.clf() 
        ECG.plot_bsqisheatmap(hamilton_rodrigues, hamilton_nabian, rodrigues_nabian)
        plt.savefig(args.file + '/bsqi/' + str(flist[i]) + 'heatmap-bsqis.png') 
        plt.clf()  
        # RESULTS & STATS for each data file of a patient -----------------------------------------------------        
        print('File path: ' + args.file + "/" + str(flist[i]))
        print('#detected R-peaks by Hamilton: ' + str(len(hamilton[1])))
        print('#detected R-peaks by Rodrigues: ' + str(len(rodrigues[1])))
        print('#detected R-peaks by Nabian: ' + str(len(nabian[1])))
        print('Length and #samples: ' + str(length) + "s , " + str(len(samples)) + " samples" )
        print('bSQI - Hamilton & Rodrigues: ' + str(hamilton_rodrigues))
        print("Average bSQI: " + str(np.round(np.average(hamilton_rodrigues),3)))
        print('Very good samples (0.8-1.0): ' + str(hrgut))
        print('Good samples (0.5-0.8): ' + str(hrmittel))
        print('Bad samples (0-0.5): ' + str(hrschlecht))
        print('Both methods detected gaps in ' + str(len(hrgaps)) + ' (' + str(hrgapscountpct)  + '%) ' + ' out of ' + str(len(samples)) + ' 10s-samples (bSQI could be corrupted)')
        print('bSQI - Hamilton & Nabian: ' + str(hamilton_nabian))
        print("Average SQI: " + str(np.round(np.average(hamilton_nabian),3)))
        print('Very good samples (0.8-1.0): ' + str(hngut))
        print('Good samples (0.5-0.8): ' + str(hnmittel))
        print('Bad samples (0-0.5): ' + str(hnschlecht))
        print('Both methods detected gaps in ' + str(len(hngaps)) + ' (' + str(hngapscountpct) + '%) ' + ' out of ' + str(len(samples)) + ' 10s-samples (bSQI could be corrupted)')
        print('bSQI - Rodrigues & Nabian: ' + str(rodrigues_nabian))
        print("Average SQI: " + str(np.round(np.average(rodrigues_nabian),3)))
        print('Very good samples (0.8-1.0): ' + str(rngut))
        print('Good samples (0.5-0.8): ' + str(rnmittel))
        print('Bad samples (0-0.5): ' + str(rnschlecht))
        print('Both methods detected gaps in ' + str(len(rngaps)) + ' (' + str(rngapscountpct) + '%) ' + ' out of ' + str(len(samples)) + ' 10s-samples (bSQI could be corrupted)')
# OVERVIEW PLOT all data files of one patient ------------------------ 
    fig1, ax = plt.subplots(datafiles, 1, figsize = (35,40))
    plt.rcParams['font.size'] = '25'  
    for i in range(datafiles):
        data = args.file + "/" + str(flist[i])
        ecg = ECG(data) 
        windowsamples = len(ecg.datasamples())
        hamilton = ecg.hamilton()
        rodrigues = ecg.rodrigues()
        nabian = ecg.nabian()
        hamilton_rodrigues = ecg.bsqi(hamilton,rodrigues)
        hamilton_nabian = ecg.bsqi(hamilton,nabian)
        rodrigues_nabian = ecg.bsqi(rodrigues,nabian)
        bsqis = []
        bsqis.append(sorted(hamilton_rodrigues))
        bsqis.append(sorted(hamilton_nabian))
        bsqis.append(sorted(rodrigues_nabian))
        ylabels = ['a)', 'b)', 'c)'] 
        if windowsamples <= 10:
            hm = sns.heatmap(bsqis,cbar_kws={'ticks': [0.0, 0.25, 0.5, 0.75, 1.0]}, vmin=0, vmax=1, xticklabels = 2, yticklabels=ylabels,ax=ax[i])
        else:
            hm = sns.heatmap(bsqis,cbar_kws={'ticks': [0.0, 0.25, 0.5, 0.75, 1.0]}, vmin=0, vmax=1, xticklabels = 20, yticklabels=ylabels,ax=ax[i])
        hm.set_xticklabels(labels=hm.get_xticklabels(),fontsize = 40)
        hm.set_yticklabels(labels=hm.get_yticklabels(), va='center', rotation = 0, fontsize = 45)
        plt.yticks(rotation=0) 
        hm.set_title(str(flist[i]))
        plt.draw()
    fig1.tight_layout()
    plt.savefig(args.file + '/bsqi/heatmap-overview.png')
    plt.clf()    
# OVERVIEW RESULTS & STATS all data files of one patient ------------------------ 
    print('--------------------------------------------------------')
    print('Total number of data files: ' + str(datafiles))
    print('Total length of all data files: ' + str(length_total))
    print('Total number of samples: ' + str(samples_total))
    print('Total number of detected R-peaks by Hamilton: ' + str(hamilton_total))
    print('Total number of detected R-peaks by Rodrigues: ' + str(rodrigues_total))
    print('Total number of detected R-peaks by Nabian: ' + str(nabian_total))
    avg_total = [hamiltonrodrigues_avgtotal,hamiltonnabian_avgtotal,rodriguesnabian_avgtotal]
    print('Total average bSQI: ' + str(np.average(avg_total)))
    print('H&R average bSQI: ' + str(np.average(hamiltonrodrigues_avgtotal)))
    print('H&R all samples Very good: ' + str(hrgut_total) + ' (' + str(np.round((100/samples_total) * hrgut_total,2)) + '%)')
    print('H&R all samples Good: ' + str(hrmittel_total) + ' (' + str(np.round((100/samples_total) * hrmittel_total,2)) + '%)')
    print('H&R all samples Bad: ' + str(hrschlecht_total) + ' (' + str(np.round((100/samples_total) * hrschlecht_total,2)) + '%)')
    print('H&R all samples Corrupted: ' + str(hrcorrupted_total) + ' (' + str(np.round((100/samples_total) * hrcorrupted_total,2)) + '%)')
    print('H&N average bSQI: ' + str(np.average(hamiltonnabian_avgtotal)))
    print('H&N all samples Very good: ' + str(hngut_total) + ' (' + str(np.round((100/samples_total) * hngut_total,2)) + '%)')
    print('H&N all samples Good: ' + str(hnmittel_total) + ' (' + str(np.round((100/samples_total) * hnmittel_total,2)) + '%)')
    print('H&N all samples Bad: ' + str(hnschlecht_total) + ' (' + str(np.round((100/samples_total) * hnschlecht_total,2)) + '%)')
    print('H&N all samples Corrupted: ' + str(hncorrupted_total) + ' (' + str(np.round((100/samples_total) * hncorrupted_total,2)) + '%)')
    print('R&N average bSQI: ' + str(np.average(rodriguesnabian_avgtotal)))
    print('R&N all samples Very good: ' + str(rngut_total) + ' (' + str(np.round((100/samples_total) * rngut_total,2)) + '%)')
    print('R&N all samples Good: ' + str(rnmittel_total) + ' (' + str(np.round((100/samples_total) * rnmittel_total,2)) + '%)')
    print('R&N all samples Bad: ' + str(rnschlecht_total) + ' (' + str(np.round((100/samples_total) * rnschlecht_total,2)) + '%)')
    print('R&N all samples Corrupted: ' + str(rncorrupted_total) + ' (' + str(np.round((100/samples_total) * rncorrupted_total,2)) + '%)')
