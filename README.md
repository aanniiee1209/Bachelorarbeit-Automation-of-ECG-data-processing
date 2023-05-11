# Automation-of-ECG-data-processing

## Imports
* Matplotlib
* Numpy
* Pandas
* OS
* fnmatch
* argparse
* Neurokit2
* datetime, timedelta
* Seaborn

## Verwendung
### bSQI-plots-stats.py
* Erstellt automatisch die Plots (Scatterplots, Histogramm, Heatmaps) für alle Messungen eines Patienten und weitere Informationen zu einzelnen und allen Messungen
  - Plots

  * für einzelne Messungen
    - File path
    - #detected R-peaks by Hamilton
    - #detected R-peaks by Rodrigues
    - #detected R-peaks by Nabian
    - Length and #samples
    - bSQI - Hamilton & Rodrigues
    - bSQI - Hamilton & Nabian
    - bSQI - Rodrigues & Nabian
    - Für jede der drei Kombination: Average bSQI, Window classification (Very good, Good, Bad), #detected gaps

  * für alle Messungen eines Patienten
    - Total number of data files
    - Total length of all data files
    - Total number of samples
    - Total number of detected R-peaks by Hamilton
    - Total number of detected R-peaks by Rodrigues
    - Total number of detected R-peaks by Nabian
    - Total average bSQI
    - für jede der drei Kombinationen: Average bSQI, Window classification (Very good, Good, Bad), #detected gaps

* Befehl für Terminal im Ordner mit .py file
```
python bSQI-plots-stats.py -f <Pfad zum Ordner mit allen Messungen eines Patienten>

```
* Beispiel
```
python bSQI-plots-stats.py -f "/Users/.../data/ecg only/pt004"

```

* Beispiel Output
```
File path: /Users/.../data/ecg only/test/17.csv
#detected R-peaks by Hamilton: 4267
#detected R-peaks by Rodrigues: 4281
#detected R-peaks by Nabian: 4395
Length and #samples: 4411.0799375s , 441 samples
bSQI - Hamilton & Rodrigues: [0.8333 1.     1.   ...]
Average bSQI: 0.979
Very good samples (0.8-1.0): 434
Good samples (0.5-0.8): 7
Bad samples (0-0.5): 0
Both methods detected gaps in 52 (11.791%)  out of 441 10s-samples (bSQI could be corrupted)
bSQI - Hamilton & Nabian: [0.8    1.     0.9091 ...]
Average SQI: 0.945
Very good samples (0.8-1.0): 416
Good samples (0.5-0.8): 24
Bad samples (0-0.5): 1
Both methods detected gaps in 29 (6.576%)  out of 441 10s-samples (bSQI could be corrupted)
bSQI - Rodrigues & Nabian: [0.6667 1.     0.9091 ...]
Average SQI: 0.962
Very good samples (0.8-1.0): 427
Good samples (0.5-0.8): 13
Bad samples (0-0.5): 1
Both methods detected gaps in 25 (5.669%)  out of 441 10s-samples (bSQI could be corrupted)

File path: /Users/.../data/ecg only/test/15.csv
#detected R-peaks by Hamilton: 5538
#detected R-peaks by Rodrigues: 5538
#detected R-peaks by Nabian: 5577
Length and #samples: 4497.62s , 449 samples
bSQI - Hamilton & Rodrigues: [0.5714 0.8125 1. ...]
Average bSQI: 0.987
Very good samples (0.8-1.0): 447
Good samples (0.5-0.8): 2
Bad samples (0-0.5): 0
Both methods detected gaps in 27 (6.013%)  out of 449 10s-samples (bSQI could be corrupted)
bSQI - Hamilton & Nabian: [0.4286 0.8125 1. ...]
Average SQI: 0.971
Very good samples (0.8-1.0): 440
Good samples (0.5-0.8): 8
Bad samples (0-0.5): 1
Both methods detected gaps in 18 (4.009%)  out of 449 10s-samples (bSQI could be corrupted)
bSQI - Rodrigues & Nabian: [0.8    1.     1. ...]
Average SQI: 0.984
Very good samples (0.8-1.0): 445
Good samples (0.5-0.8): 4
Bad samples (0-0.5): 0
Both methods detected gaps in 17 (3.786%)  out of 449 10s-samples (bSQI could be corrupted)
--------------------------------------------------------
Total number of data files: 2
Total length of all data files: 8908.6999375
Total number of samples: 890
Total number of detected R-peaks by Hamilton: 9805
Total number of detected R-peaks by Rodrigues: 9819
Total number of detected R-peaks by Nabian: 9972
Total average bSQI: 0.9713333333333334
H&R average bSQI: 0.983
H&R all samples Very good: 881 (98.99%)
H&R all samples Good: 9 (1.01%)
H&R all samples Bad: 0 (0.0%)
H&R all samples Corrupted: 79 (8.88%)
H&N average bSQI: 0.958
H&N all samples Very good: 856 (96.18%)
H&N all samples Good: 32 (3.6%)
H&N all samples Bad: 2 (0.22%)
H&N all samples Corrupted: 47 (5.28%)
R&N average bSQI: 0.973
R&N all samples Very good: 872 (97.98%)
R&N all samples Good: 17 (1.91%)
R&N all samples Bad: 1 (0.11%)
R&N all samples Corrupted: 42 (4.72%)

```

### bSQI-ecg.py
* Öffnet EKG mit 10 s-Fenster Einteilung, R-peak detections der zwei verwendeten Methoden und den bSQI Werten der Fenster

* Befehl für Terminal im Ordner mit .py file
```
python bSQI-ecg.py -f <Pfad zu einer Messungen eines Patienten> -m1 <Methode 1> -m2 <Methode 2>

```
* Beispiel
```
python bSQI-ecg.py -f "/Users/.../data/ecg only/pt007/movesense_ecg_20220920130534933.csv" -m1 "hamilton" -m2 "rodrigues"

```

 - Wenn nur die R-peak detection einer Methode gezeigt werden sollen, dann für -m1 und -m2 die gleiche Methode eintragen (in diesem Fall ist der bSQI für alle Fenster = 1)
* Beispiel
```
python bSQI-ecg.py -f "/Users/.../data/ecg only/pt007/movesense_ecg_20220920130534933.csv" -m1 "hamilton" -m2 "hamilton"

```

