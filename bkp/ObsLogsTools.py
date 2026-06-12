def getObsLogForCurrentScan(doReturn=False):
   """Method to search for ObsLog info of current scan. 
   For this to work the collection of obslog files of the current account 
   has to be saved in the folder:~/obslogs/
   """

   scan=data.ScanParam.ScanNum
   day=data.ScanParam.DateObs[:10]
   find_log=False
   files= os.listdir('obslogs')
   for infile in files:
      if day in infile:
         f=open('obslogs/'+infile,'r')
         lines=f.readlines()
         f.close()
         find_log=True
         break
   if not find_log:
      print('No obslog for %s'%day)
      return
   keys=[]
   index=0
   find_scan=False
   message=''
   for index in range(len(lines)):
      line=lines[index]
      if line[0:4]=='<th>':
         keys.append(line[4:-6])
         index+=1
      if line[0:4]=='<td>' and line[4:-6]==str(scan):
         find_scan=True
         index+=1
         for key in keys:
            line=lines[index]
            message += (key + ' , ' +line[4:-6] + '\n')
            index+=1
         break
   if not find_scan:
      print('No entry for scan %i in %s'%(scan,day))
   if doReturn:
      return message
   else:
      print(meassage)
      return

def ScansOverview(Scans=None):
   """Method to search for ObsLog info of current account.
   It display an overview with basic information. For more detais use:
   getObsLogForCurrentScan().   
   For this method to work the collection of obslog files of the current account 
   has to be saved in the folder:~/obslogs/
   """ 
   files= os.listdir('obslogs')
   for infile in files:
      __getScansOverview('obslogs/'+infile,Scans=Scans)
   return
   
def __getScansOverview(infile,Scans=None):
  
   f=open(infile,'r')
   lines=f.readlines()
   index=0
   start=False
   keys=[]
   for index in range(len(lines)):
      line=lines[index]
      if line[0:4]=='<th>':
         keys.append(line[4:-6])
         index+=1
      elif line[0:4]=='<tr>':
         start=True
         index+=1    
      elif line[0:5]=='</tr>':
         index+=1
      else:
         index+=1
      if start:   
         message=''   
         scan=0 
         for key in keys:
            line=lines[index]
            index+=1 
            if key=='Scan' :
               scan=int(line[4:-6])
               message+=(line[4:-6] + ' , ')
            if key=='Source':
               message+=(line[4:-6].ljust(12) + ' , ')               
            if key== 'Scan status':
               message+=(line[4:-6] + ' , ')
            if key=='Scan type':
               message+=(line[4:-6].ljust(6) + ' , ')
            if key=='Comment':
               message+=(line[4:-6].ljust(20) )               
         start=False
         if Scans:
            if 'OK' in message and scan in Scans:
               print(message)
         else:
            if 'OK' in message:
               print(message)

   return

