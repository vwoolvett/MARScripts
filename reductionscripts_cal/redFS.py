f = open('CalFiles/allFS.dat')
lines = f.readlines()
f.close()

scans=[]
for line in lines:
    scans.append(np.int(line.split("_")[2]))
   
for scan in scans:
    try:
        reduceFsweep(scan, overwrite=True, doPlot=False)
    except:
        continue 
