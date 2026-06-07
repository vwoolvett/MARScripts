import copy as copy

scans = []
bad = []
for s in bad:
    if s in scans:
        scans.remove(s)


myname = 'LFA-G345'
fe = 'LFA'
writeSummary=True
clip = 3.

##map size  in absolute EQ or GAL coordinates in deg :
system='EQ'
xsize=[256.7,252.0]
ysize=[-42.0,-38.0]
    
doPlot=False
if not doPlot:
    noPlot=True
else:
    noPlot=False


for iter in range(1,3):
    if iter == 1:
        mymodel=None
        subtract = False
    else:
        mymodel=str(myname)+"-coadded-flux-iter"+str(iter-1)+".data"
        m=restoreFile(mymodel)
        if iter == 2:
            subtract = False
            mymodel=createSourceModel(m,highcut=5.5,lowcut=2.5,sm=8.,mtype='snr',clip=3)
            
        if iter == 3:
            subtract = True
            mymodel=createSourceModel(m,highcut=5.5,lowcut=2.5,sm=0.,mtype='flux',clip=3)
            
    
    ms=None

    for i,scan in enumerate(scans):
        scanname = "ReducedFiles/"+str(myname)+"-"+str(scan)+"-iter"+str(iter)+".data"
        print scanname
        globlist=glob(scanname)

        m = None
        if len(globlist) ==  0:
            redweak(scan,fe='LFA',size=-1,model=mymodel,subtract=subtract,doPlot=doPlot,extremeFilter=True,writeSummary=writeSummary)
            #if scan == 22919: #flagging example to flag a certain time range in a map (seconds from the beining of the scan) 
            #    flagMJD(above=1430,below=1600,flag=2)
            mapping(oversamp=4,system=system,sizeX=xsize,sizeY=ysize,limitsZ=[-0.8,1.5],noPlot=noPlot)
            data.Map.dumpMap(scanname)
            m=restoreFile(scanname)
            m.smoothBy(8./3600.)
            if writeSummary:
                rmsMap = copy.deepcopy(m)
                rmsMap.Data =  (rmsMap.Data*0.0+1.0)/np.sqrt(rmsMap.Weight)
                minnoise=np.nanmin(rmsMap.Data)
                mask=np.where(rmsMap.Data > 5*minnoise)
                rmsMap.Data[mask] = np.NaN
                pixelsize=np.abs(rmsMap.WCS['CDELT2'])
                nrpix = np.sum(~np.isnan(rmsMap.Data))
                area = nrpix*pixelsize**2
                noise=np.nanmedian(rmsMap.Data)
                outname="%s-%s-%i_summary.txt"%(fe,data.ScanParam.Object,data.ScanParam.ScanNum)
                f=open(outname,'r')
                lines=f.readlines()
                f.close()
                myline=lines[0].replace("\n","")
                myline=myline+",{:.3f},{:.4f}\n".format(area,noise)
                f=open(outname,'w')
                f.write(myline)
                f.close()
            
        else:
            m=restoreFile(scanname)
            m.smoothBy(8./3600.)
            

        if ms and m:
            ms = mapsumfast([ms,m])
        elif not ms:
            ms = copy.deepcopy(m)
    
        
    
        rmsMap = copy.deepcopy(ms)
        snrMap = copy.deepcopy(ms)

        snrMap.Data *= np.sqrt(snrMap.Weight)
        a = snrMap.computeRms()
        scale=snrMap.RmsBeam
        snrMap.Data /= np.array(scale,'f')

        rmsMap.Data =  (rmsMap.Data*0.0+1.0)/np.sqrt(rmsMap.Weight)
        rmsMap.Data *= np.array(scale,'f')
    
        if doPlot:
            snrMap.display(aspect=1,limitsZ=[-4,12])

    minnoise=np.nanmedian(rmsMap.Data)
    meannoise=np.nanmedian(rmsMap.Data)
    
    if clip > 0:
        mask=np.where(rmsMap.Data > 5*minnoise)
        ms.Data[mask] = np.NaN
        snrMap.Data[mask] = np.NaN
        rmsMap.Data[mask] = np.NaN

    snrMap.display(aspect=0,limitsZ=[-4,12])
    rmsMap.display(aspect=0,limitsZ=[0,2*meannoise],doContour=1,levels=[meannoise],overplot=1)

    outname=str(myname)+"-coadded-flux-iter"+str(iter)+".data"
    ms.dumpMap(outname)


    print("###################iteration %i ##########################"%(iter))
    print("minimum noise: %5.1f mJy/b, mean noise: %5.1f mJy/b"%(1000*minnoise,1000*meannoise))
    print("#####################################################")



    outname=str(myname)+"-coadded-iter"+str(iter)+".fits"
    writeFits2(ms,outfile=outname,overwrite=1)

