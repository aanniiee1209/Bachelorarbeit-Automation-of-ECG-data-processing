# Automation-of-ECG-data-processing

## Needed imports
* Matplotlib
* Numpy
* Pandas
* OS
* fnmatch
* argparse
* Neurokit2
* datetime, timedelta
* Seaborn

## Usage
### bSQI-plots-stats.py
* Erstellt die Plots (Scatterplots, Histogramm, Heatmap) für alle Messungen eines Patienten und weitere Informationen zu einzelnen Files und allen Files
 - Plots

 * für einzelne Files
  - File path
  - #detected R-peaks by Hamilton
  - #detected R-peaks by Rodrigues
  - #detected R-peaks by Nabian
  - Length and #samples
  - bSQI - Hamilton & Rodrigues
  - bSQI - Hamilton & Nabian
  - bSQI - Rodrigues & Nabian
  - Für jede Kombination: Average bSQI, Window classification (Very good, Good, Bad), #detected gaps by both methods

 * für alle Messungen eines Patienten
  - Total number of data files
  - Total length of all data files
  - Total number of samples
  - Total number of detected R-peaks by Hamilton
  - Total number of detected R-peaks by Rodrigues
  - Total number of detected R-peaks by Nabian
  - Total average bSQI
  - für jede Kombination: Average bSQI, Window classification (Very good, Good, Bad), #detected gaps by both methods

* Example Output
```
File path: /Users/ngocdonganhvo/Library/CloudStorage/OneDrive-bwedu/Bachelor Unizeug/Bachelorarbeit/data/ecg only/test/15.csv
#detected R-peaks by Hamilton: 5538
#detected R-peaks by Rodrigues: 5538
#detected R-peaks by Nabian: 5577
Length and #samples: 4497.62s , 449 samples
bSQI - Hamilton & Rodrigues: \[...]
Average bSQI: 0.987
Very good samples (0.8-1.0): 447
Good samples (0.5-0.8): 2
Bad samples (0-0.5): 0
Both methods detected gaps in 27 (6.013%)  out of 449 10s-samples (bSQI could be corrupted)
bSQI - Hamilton & Nabian: \[...]
Average SQI: 0.971
Very good samples (0.8-1.0): 440
Good samples (0.5-0.8): 8
Bad samples (0-0.5): 1
Both methods detected gaps in 18 (4.009%)  out of 449 10s-samples (bSQI could be corrupted)
bSQI - Rodrigues & Nabian: \[...]
Average SQI: 0.984
Very good samples (0.8-1.0): 445
Good samples (0.5-0.8): 4
Bad samples (0-0.5): 0 
Both methods detected gaps in 17 (3.786%)  out of 449 10s-samples (bSQI could be corrupted)

--------------------------------------------------------
Total number of data files: 3
Total length of all data files: 11977.687937499999
Total number of samples: 1197
Total number of detected R-peaks by Hamilton: 13788
Total number of detected R-peaks by Rodrigues: 13801
Total number of detected R-peaks by Nabian: 13956
Total average bSQI: 0.9755555555555555
H&R average bSQI: 0.9843333333333333
H&R all samples Very good: 1187 (99.16%)
H&R all samples Good: 9 (0.75%)
H&R all samples Bad: 1 (0.08%)
H&R all samples Corrupted: 79 (6.6%)
H&N average bSQI: 0.965
H&N all samples Very good: 1161 (96.99%)
H&N all samples Good: 33 (2.76%)
H&N all samples Bad: 3 (0.25%)
H&N all samples Corrupted: 47 (3.93%)
R&N average bSQI: 0.9773333333333333
R&N all samples Very good: 1177 (98.33%)
R&N all samples Good: 18 (1.5%)
R&N all samples Bad: 2 (0.17%)
R&N all samples Corrupted: 42 (3.51%)
```

* Befehl für Terminal im Ordner mit .py file
```
python bSQI-plots-stats.py -f <Pfad zum Ordner mit allen Messungen eines Patienten>

```
* Example
```
python bSQI-plots-stats.py -f "/Users/.../data/ecg only/pt004"

```

### bSQI-ecg.py
* Befehl für Terminal im Ordner mit .py file
```
python bSQI-ecg.py -f <Pfad zu einer Messungen eines Patienten> -m1 <Methode 1> -m2 <Methode 2>

```
* Example
```
python bSQI-ecg.py -f "/Users/.../data/ecg only/pt007/movesense_ecg_20220920130534933.csv" -m1 "hamilton" -m2 "rodrigues"

```
 - Wenn nur die R-peak detection einer Methode gezeigt werden sollen, dann für -m1 und -m2 die gleiche Methode eintragen (in diesem Fall ist der bSQI für alle Fenster = 1)
* Example
```
python bSQI-ecg.py -f "/Users/.../data/ecg only/pt007/movesense_ecg_20220920130534933.csv" -m1 "hamilton" -m2 "hamilton"

```

