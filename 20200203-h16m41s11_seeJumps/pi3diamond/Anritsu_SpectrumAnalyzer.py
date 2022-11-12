from __future__ import print_function, absolute_import, division


import visa
rm = visa.ResourceManager()
anritsu = rm.open_resource('GPIB0::11::INSTR', timeout=20000)

def SetCenterFreq(freq):
    try:
        f=float(freq)
        anritsu.write("CF "+str(f)+"HZ")
    except:
        pass
    
def SetSpan(freq):
    try:
        f=float(freq)
        anritsu.write("SP "+str(f)+"HZ")
    except:
        pass
    
def doSweep():
    anritsu.write("TS")
    
def getPeak():
    doSweep()
    anritsu.write("MKPK")
    freq=float(anritsu.ask("MKF?"))
    amp=float(anritsu.ask("MKL?"))
    return freq, amp
    
def getDip():
    doSweep()
    anritsu.write("MKMIN")
    freq=float(anritsu.ask("MKF?"))
    amp=float(anritsu.ask("MKL?"))
    return freq, amp