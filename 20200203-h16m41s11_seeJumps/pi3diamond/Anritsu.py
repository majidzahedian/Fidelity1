"""
This file is part of Diamond. Diamond is a confocal scanner written
in python / Qt4. It combines an intuitive gui with flexible
hardware abstraction classes.

Diamond is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Diamond is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with diamond. If not, see <http://www.gnu.org/licenses/>.

Copyright (C) 2009 Helmut Rathgen <helmut.rathgen@gmail.com>
"""

from __future__ import print_function, absolute_import, division

import visa
import time
# get device string from pi3d

class Anritsu(object):

    def __init__(self):
        rm = visa.ResourceManager(r"C:\Windows\System32\visa64.dll")
        self.anritsu = rm.open_resource('TCPIP::192.168.0.2::inst0::INSTR', timeout=20000)
    

    def On(self):
        self.anritsu.write('OUTP:STAT ON')
        self.anritsu.write('*WAI')
        return self.anritsu.query('OUTP:STAT?')

    def Off(self):
        self.anritsu.write('OUTP:STAT OFF')
        t0 = time.time()
        try:
            return self.anritsu.query('OUTP:STAT?')
        except Exception as e:
            print('It took {} seconds'.format(time.time() - t0))
            raise e

    def Power(self, power=None):
        if power != None:
            self.anritsu.write(':POW %f' % power)
        return float(self.anritsu.query(':POW?'))

    def Freq(self, f=None):
        if f != None:
            self.anritsu.write(':FREQ %f' % f)
        return float(self.anritsu.query(':FREQ?'))

    def CW(self, f=None, power=None):
        self.anritsu.write(':FREQ:MODE CW')
        if f != None:
            self.anritsu.write(':FREQ %f' % f)
        if power != None:
            self.anritsu.write(':POW %f' % power)

    def List(self, freq, power, start_pos=0):
        self.anritsu.write(':POW %f' % power)
        self.anritsu.write(':LIST:TYPE FREQ')
        self.anritsu.write(':LIST:IND 0')
        s = ''
        for f in freq[:-1]:
            s += ' %f,' % f
        s += ' %f' % freq[-1]
        self.anritsu.write(':LIST:FREQ' + s)
        self.anritsu.write(':LIST:STAR 0')
        self.anritsu.write(':LIST:STOP %i' % (len(freq)-1) )
        self.anritsu.write(':LIST:MODE MAN')
        self.anritsu.write(':LIST:IND %s'%start_pos)
        self.anritsu.write('*WAI')

    def List_fp(self, freq, power):
        #self.anritsu.write(':POW %f' % power)
        self.anritsu.write(':LIST:TYPE FLEV')
        self.anritsu.write(':LIST:IND 0')
        s = ''
        for f in freq[:-1]:
            s += ' %f,' % f
        s += ' %f' % freq[-1]
        self.anritsu.write(':LIST:FREQ' + s)
        self.anritsu.write(':LIST:IND 0')
        s = ''
        for p in power[:-1]:
            s += ' %f,' % p
        s += ' %f' % power[-1]
        self.anritsu.write(':LIST:POW' + s)
        self.anritsu.write(':LIST:STAR 0')
        self.anritsu.write(':LIST:STOP %i' % (len(freq)-1) )
        self.anritsu.write(':LIST:MODE MAN')
        self.anritsu.write('*WAI')

    def Trigger(self, source, pol):
        self.anritsu.write(':TRIG:SOUR '+source)
        self.anritsu.write(':TRIG:SLOP '+pol)
        self.anritsu.write('*WAI')

    def ResetListPos(self):
        self.anritsu.write(':LIST:IND 0')
        self.anritsu.write('*WAI')

    def GetListPos(self):
        return int(self.anritsu.query(':LIST:IND?'))

    def ListOn(self):
        self.anritsu.write(':FREQ:MODE LIST')
        self.anritsu.write(':OUTP ON')
        self.anritsu.write('*WAI')

    def Modulation(self, flag=None):
        if flag is not None:
            if flag:
                self.anritsu.write('PULM:SOUR EXT')
                self.anritsu.write('PULM:STAT ON')
            else:
                self.anritsu.write('PULM:STAT OFF')
        return bool(int(self.anritsu.query('PULM:STAT?')))
