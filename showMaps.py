#indir('/apex-archive/rawdata/M-0117.F-9550A-2026/')
import copy as copy


scans = [27974,27975,27979,27980,27990,28212,28213,28217,28218,28231,28232,28235,28493,28494,28498,28499,28516,28517,28775,28776]
bad=[27979,28217,28498]

check=[28213]

iter=2
myname = 'LFA-OG259-EQ'
fe = 'LFA'

for i,scan in enumerate(scans):
    scanname = "ReducedFiles/"+str(myname)+"-"+str(scan)+"-iter"+str(iter)+".data"
    print scanname
    globlist=glob(scanname)

      
    m=restoreFile(scanname)
    m.smoothBy(8./3600.)
    m.display(aspect=1,limitsZ=[-0.2,0.5])    
    raw_input()
    
