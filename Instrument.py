import visa
import numpy as np
import time

rm = visa.ResourceManager()

#TODO setter do Attenuator não consegue adicionar db de 100


# Wrap phase______________
def unwrap_d(Phase):
    """
    Unwrap Phase in Degree"""
    return np.unwrap(Phase * np.pi / 180)


def wrap(Phase):
    """ Wrap back the phase 
    """
    return (Phase + np.pi) % (2 * np.pi) - np.pi


def wrap1(Phase):
    """ Wrap back the phase using arctan2
    """
    return arctan2(sin(Phase), cos(Phase))

    


def dBmtoV(Amp):
    """Convert dBm to Volts
    Input Amplitude in dBm
    """
    c = 10 * np.log10(20)
    return np.power(10, (Amp - c) / 20)


def VtodBm(V):
    """Convert Volts to dBm
    Input Voltage in Volts
    """
    c = 10 * np.log10(20)
    return 20 * np.log10(V) + c


# Convert dBm to W

def dBmtoW(Amp):
    """Convert dBm to Watts
    Imput Amp in dBm
    """
    return np.power(10, (Amp - 30) / 10)

def asc_to_array_int(data):
    data = data.split('\n')[0].split(',')
    return [int(float(i)) for i in data]

def asc_to_array_float(data):
    data = data.split('\n')[0].split(',')
    return [float(i) for i in data]

def fit_peak(X,Y,n):
    from scipy.optimize import curve_fit
    def func(x, a, b, c):
        return a*(x-b)**2+c

    peak = np.where(Y == np.min(Y))[0][0]
    begin = peak-n
    end = peak +n
    n_interpolate = 41
    xnew = np.linspace(X[begin], X[end-1], num=n_interpolate, endpoint=True)
    popt, pcov = curve_fit(func, X[begin:end], Y[begin:end],p0=[1,X[peak],Y[peak]])
    ynew= func(xnew,*popt)
    
    return xnew,ynew

def get_curve(na):
    na.data_array_format = 'PLOGarithmic'
    return na.data_array(),na.frequency_array()

def parse_data(data,freq):
    Y = [float(i) for i in data.split('\n')[0].split(',')]
    phase = np.array(Y[1::2])
    del Y[1::2]
    mag = np.array(Y)
    X = np.array(asc_to_array_int(freq))

    return X,mag,phase




class Instrument:
    def __init__(self,address,alias):
        resource = rm.open_resource(address)
        self._inst = resource
        self._idn = self._inst.query("*IDN?")[:-1]
        self._alias = alias
        
    def write(self,command):
        return self._inst.write(command)
    
    def write_raw(self,commad):
        return self._inst.write_raw(commad)
        
    def query(self,command):
        return self._inst.query(command)
    
    def read_raw(self):
        return self._inst.read_raw()
    
    def open(self):
        return self._inst.open()
    
    def close(self):
        return self._inst.close()
    
class AWG(Instrument):
    def __init__(self, address, alias):
        super().__init__(address,alias)
        self.stop()
        self.write(":INST:DACM SING")
        self.write(":INST:MEM:EXT:RDIV DIV1")
        
    def _convertToByte(self,data, A):

        np.clip(data, -A, A, data)

        size = 256 * (1 + divmod(len(data) - 1, 256)[0])
        
        y = np.zeros(size, dtype=np.int8)
        y[:len(data)] = np.array(127 * data / A, dtype=np.int8)
        return y
    
    def loadWaveform(self,data):
        self.convertedData = self._convertToByte(data,self.Vamp)
        
    def clearWaveform(self,ch,seq):
        self.write(":TRAC{}:DEL {}".format(ch,seq));
        
    def sendWaveform(self,ch,seq):
        data = self.convertedData
        n_elem = len(data)
        
        self.write(':TRAC{}:DEF {},{},0'.format(ch,seq, n_elem))
        # create binary data as bytes with header
        start, length = 0, len(data)
        sLen = b'%d' % length
        sHead = b'#%d%s' % (len(sLen), sLen)
        # send to AWG
        sCmd = b':TRAC%d:DATA %d,%d,' % (ch, seq, start)
        self.write_raw(sCmd + sHead + data[start:start+length].tobytes())
        
    def start(self):
        self.write('INIT:IMM')
        
    def stop(self):
        self.write(':ABOR')
        
    def enableCh(self,ch):
        self.write(":OUTP%d ON" % ch)
        
    def disableCh(self,ch):
        self.write(":OUTP%d OFF" % ch)
        
    def setTriggerToContinuous(self):
        self.write(":INIT:CONT ON")
        
    def setTriggerToGated(self):
        self.write(":INIT:GATE ON")
        
    def setTriggerToTriggered(self):
        self.write(":INIT:CONT OFF")
        self.write(":INIT:GATE OFF")
        
    def forceTrigger(self):
        awg.write(":TRIG:BEG")
    
class Attenuator(Instrument):
    def __init__(self,address,alias):
        super().__init__(address,alias)
        self._attenuation = self.attenuation
    
    @property
    def attenuation(self):
        '''This sets the attenuation level. '''
        
        d = int(self.query("ATTenuator:BANK1:Y?").split("\n")[0])
        u = int(self.query("ATTenuator:BANK1:X?").split("\n")[0])
        self._attenuation = d+u
        return self._attenuation
    
    
    @attenuation.setter
    def attenuation(self,att):
        '''This gets the attenuation level. '''
        
        d = int(att/10)
        u = int(att%10)
        self.write('ATT:BANK1:X {}'.format(u))
        self.write('ATT:BANK1:Y {}'.format(d*10))
        self._attenuation = d+u
    
class NetworkAnalyzer(Instrument):
    def __init__(self,address,alias):
        super().__init__(address,alias)
        self._average_points = self.average_points
        self._center_freq = self.center_frequency 
        self._start_freq = self.start_frequency 
        self._stop_freq = self.stop_frequency 
        self._span_freq = self.span_frequency
        self._sweep_points = self.sweep_points
        self._if_bandwidth = self.if_bandwidth_frequency
        self._data_format = self.data_format
        
    def autoscale(self):
        '''This command executes the auto scale function. The Auto Scale function automatically adjusts the value of the reference division line and the scale per division to display the trace appropriately.'''
        self.write(':DISPlay:WINDow1:TRACe1:Y:SCALe:AUTO')
        
    @property
    def sweep_time(self):
        '''This command gets the sweep time.
        
            Parameter: Value
            Description: Sweep time
            Data Type: Numeric type (Real)
            Range: 0 ~ 2M
            Preset Value: 24.7228364697 ms
            Unit: s
            Resolution: -
        '''
        return float(self.query(":SENSe1:SWEep:TIME:DATA?").split('\n')[0])
        
    @property
    def smooth(self):
        '''This command turns ON/OFF the smoothing, for the active trace of the selected channel Ch.
        
            Parameter: Selection Option
            Description: ON/OFF of the smoothing
            Data Type: Boolean type (Boolean)
            Range: ON|OFF|1|0
            Preset Value: OFF
        '''
        return self.query(":CALCulate1:SELected:SMOothing:STATe?").split("\n")[0]
    
    @smooth.setter
    def smooth(self,state):
        '''This command turns ON/OFF the smoothing, for the active trace of the selected channel Ch.
        
            Parameter: Selection Option
            Description: ON/OFF of the smoothing
            Data Type: Boolean type (Boolean)
            Range: ON|OFF|1|0
            Preset Value: OFF
        '''
        self.write(":CALCulate1:SELected:SMOothing:STATe {}".format(state))
        
        
    @property
    def power(self):
        '''This command turns ON/OFF of the stimulus signal output. Measurement cannot be made until the stimulus signal output is turned ON.
        
            Parameter: Selection Option
            Description: Status of the signal output
            Data Type: Boolean type (Boolean)
            Range: ON|OFF|1|0
            Preset Value: ON
        '''
        return self.query('OUTPut:STATe?').split('\n')[0]
    @power.setter
    def power(self,state):
        '''This command turns ON/OFF of the stimulus signal output. Measurement cannot be made until the stimulus signal output is turned ON.
        
            Parameter: Selection Option
            Description: Status of the signal output
            Data Type: Boolean type (Boolean)
            Range: ON|OFF|1|0
            Preset Value: ON
        '''
        self.write('OUTPut:STATe {}'.format(state))
        
    
    @property
    def port(self):
        '''This command gets the output port used for absolute.
        
            Parameter: Value
            Description: Source output port
            Data Type: Numeric type (Integer)
            Range: 1 to 2
            Preset Value: 1
            Note: You need to set the measurement parameter for absolute measurements with the   :CALC:PAR:DEF command.
        '''
        return int(self.query(":CALCulate1:PARameter1:SPORt?").split('\n')[0])
    
    @port.setter
    def port(self,p):
        '''This command sets the output port used for absolute.
        
            Parameter: Value
            Description: Source output port
            Data Type: Numeric type (Integer)
            Range: 1 to 2
            Preset Value: 1
            Note: You need to set the measurement parameter for absolute measurements with the   :CALC:PAR:DEF command.
        '''
        self.write(":CALCulate1:PARameter1:SPORt {}".format(p))
    
    @property
    def parameter(self):
        '''This command gets the measurement parameter. :CALC:PAR:SPOR specifies the source output port for the absolute measurements.
        
            Parameter: Selection Option
            Description: Measurement parameter
            Data Type: Character string type (String)
            Range: Select either one of the following:
                    "S<XY>"
                    Where:
                    x=1 to 2
                    Y=1 to 2
            Preset Value: "S11"
        '''
        return self.query(":CALCulate1:PARameter1:DEFine?").split('\n')[0]
    
    @parameter.setter
    def parameter(self,par):
        '''This command gets the measurement parameter. :CALC:PAR:SPOR specifies the source output port for the absolute measurements.
        
            Parameter: Selection Option
            Description: Measurement parameter
            Data Type: Character string type (String)
            Range: Select either one of the following:
                    "S<XY>"
                    Where:
                    x=1 to 2
                    Y=1 to 2
            Preset Value: "S11"
        '''
        
        self.write(":CALCulate1:PARameter1:DEFine {}".format(par))
        
    def data_array(self):
        '''This command sets/gets the formatted data array.

            The array data element varies in the data format. For more information, see Internal Data Processing.

            If valid data is not calculated because of the invalid measurement, “1.#QNB” is read out.
            
            Parameter: Value
            Description: Formatted data array

                        Where n is an integer between 1 and NOP (number of measurement points):

                            <numeric n×2-1>:  Real part of data (complex number) at the n-th measurement point.

                            <numeric n×2>: Imaginary part of data (complex number) at the n-th measurement point. Always 0 when the data format is not the Smith chart format or the polar format

                        The number of data is {NOP×2} 
            Data Type: Variant type Array (Range) 
            Note: If there is no array data of NOP×2 when setting a formatted data array, an error occurs when executed.
        '''
        
        
        data = self.query(":CALCulate1:SELected:DATA:FDATa?")
        return data
        
    def frequency_array(self):
        '''This command returns the frequency stimulus data.
        
        Parameter: Value
        Description: Indicates the array data (frequency) of NOP (number of measurement points). Where n is an integer between 1 and NOP.

                        Data(n-1): Frequency at the n-th measurement point

                    The index of the array starts from 0. 
        Data Type: Variant type (Variant)'''
        
        
        data = self.query("SENSe:FREQuency:DATA?")
        return data
    
    @property
    def data_array_format(self):
        '''This command gets the data array format.
        
            Parameter: Selection Option
            Description: Data format
            Data Type: Character string type (String)
            Range: "MLOGarithmic": Specifies the log magnitude format.
                    "PHASe": Specifies the phase format.
                    "GDELay": Specifies the group delay format.
                    "SLINear": Specifies the Smith chart format (Lin/Phase).
                    "SLOGarithmic": Specifies the Smith chart format (Log/Phase).
                    "SCOMplex": Specifies the Smith chart format (Re/Im).
                    "SMITh": Specifies the Smith chart format (R+jX).
                    "SADMittance": Specifies the Smith chart format (G+jB).
                    "PLINear": Specifies the polar format (Lin/Phase).
                    "PLOGarithmic": Specifies the polar format (Log/Phase).
                    "POLar": Specifies the polar format (Re/Im).
                    "MLINear": Specifies the linear magnitude format.
                    "SWR": Specifies the SWR format.
                    "REAL": Specifies the real format.
                    "IMAGinary": Specifies the imaginary format.
                    "UPHase": Specifies the expanded phase format.
                    "PPHase": Specifies the positive phase format. 
            Preset Value: "MLOGarithmic"
        '''
        return self.query(":CALCulate1:SELected:FORMat?").split('\n')[0]
    @data_array_format.setter
    def data_array_format(self,dformat):
        '''This command sets the data array format.
        
            Parameter: Selection Option
            Description: Data format
            Data Type: Character string type (String)
            Range: "MLOGarithmic": Specifies the log magnitude format.
                    "PHASe": Specifies the phase format.
                    "GDELay": Specifies the group delay format.
                    "SLINear": Specifies the Smith chart format (Lin/Phase).
                    "SLOGarithmic": Specifies the Smith chart format (Log/Phase).
                    "SCOMplex": Specifies the Smith chart format (Re/Im).
                    "SMITh": Specifies the Smith chart format (R+jX).
                    "SADMittance": Specifies the Smith chart format (G+jB).
                    "PLINear": Specifies the polar format (Lin/Phase).
                    "PLOGarithmic": Specifies the polar format (Log/Phase).
                    "POLar": Specifies the polar format (Re/Im).
                    "MLINear": Specifies the linear magnitude format.
                    "SWR": Specifies the SWR format.
                    "REAL": Specifies the real format.
                    "IMAGinary": Specifies the imaginary format.
                    "UPHase": Specifies the expanded phase format.
                    "PPHase": Specifies the positive phase format. 
            Preset Value: "MLOGarithmic"
        '''
        self.write(":CALCulate1:SELected:FORMat {}".format(dformat))
        
    
    @property
    def data_format(self):
        '''This command can be used to set/get the format data using the following SCPI commands:
            :CALC:DATA:FDAT
            :CALC:DATA:FMEM
            :CALC:DATA:SDAT
            :CALC:DATA:SMEM
            :CALC:DATA:XAX
            :CALC:FUNC:DATA
            :CALC:LIM:DATA
            :CALC:LIM:REP
            :CALC:LIM:REP:ALL
            :CALC:BLIM:REP
            :CALC:PLIM:DATA
            :CALC:PLIM:REP
            :CALC:RLIM:DATA
            :CALC:PLIM:REP
            :SENS:CORR:COEF
            :SENS:FREQ:DATA
            :SENS:SEGM:DATA

            Parameter: Selection Option
            Description: Data transfer format
            Data Type: Character string type (String)
            Range: ASCii: ASCII transfer format
                   REAL: IEEE 64-bit floating point binary transfer format
                   REAL32: IEEE 32-bit floating point binary transfer format
            Preset Value: ASCii
        '''
        dformat = self.query(":FORMat:DATA?").split('\n')[0]
        self._data_format = dformat
        return self._data_format
    @data_format.setter
    def data_format(self,dtype):
        '''This command can be used to sett the format data using the following SCPI commands:
            :CALC:DATA:FDAT
            :CALC:DATA:FMEM
            :CALC:DATA:SDAT
            :CALC:DATA:SMEM
            :CALC:DATA:XAX
            :CALC:FUNC:DATA
            :CALC:LIM:DATA
            :CALC:LIM:REP
            :CALC:LIM:REP:ALL
            :CALC:BLIM:REP
            :CALC:PLIM:DATA
            :CALC:PLIM:REP
            :CALC:RLIM:DATA
            :CALC:PLIM:REP
            :SENS:CORR:COEF
            :SENS:FREQ:DATA
            :SENS:SEGM:DATA

            Parameter: Selection Option
            Description: Data transfer format
            Data Type: Character string type (String)
            Range: ASCii: ASCII transfer format
                   REAL: IEEE 64-bit floating point binary transfer format
                   REAL32: IEEE 32-bit floating point binary transfer format
            Preset Value: ASCii
        '''
        self.write(":FORMat:DATA {}".format(dtype))
        
        
    def clear_averaging(self):
        '''This command resets the data count to 0, it is used for averaging. Measurement data before the execution of this object is not used for averaging.'''
        self.write(":SENSe1:AVERage:CLEar")
        
    
    @property
    def averaging(self):
        '''This command gets the averaging function.
        
            Parameter: Selection Option
            Description: ON/OFF of the averaging function
            Data Type: Boolean type (Boolean)
            Range: ON|OFF|1|0
            Preset Value: OFF
        '''
        return self.query(":SENSe1:AVERage:STATe?").split('\n')[0]
        
    @averaging.setter
    def averaging(self,state):
        '''This command sets the averaging function.
        
            Parameter: Selection Option
            Description: ON/OFF of the averaging function
            Data Type: Boolean type (Boolean)
            Range: ON|OFF|1|0
            Preset Value: OFF
        '''
        
        self.write(":SENSe1:AVERage:STATe {}".format(state))

    
    @property
    def if_bandwidth_frequency(self):
        '''This command gets the IF bandwidth.

            Parameter: Value
            Description: IF bandwidth
            Data Type: Numeric type (Real)
            Range: 10|15|20|30|40|50|70|100|150|200|300|400|500|700|1k|1.5k|2k|3k|4k|5k|7k|10k|15k|20k|30k|40k|50k|70k|100k|150k|200k|300k
            Preset Value: 70k
            Unit:Hz
        '''
        
        if_band = int(float(self.query(":SENSe1:BANDwidth:RESolution?").split("\n")[0]))
        self._if_bandwidth = if_band
        return self._if_bandwidth
    @if_bandwidth_frequency.setter
    def if_bandwidth_frequency(self,freq):
        '''This command sets the IF bandwidth.

            Parameter: Value
            Description: IF bandwidth
            Data Type: Numeric type (Real)
            Range: 10|15|20|30|40|50|70|100|150|200|300|400|500|700|1k|1.5k|2k|3k|4k|5k|7k|10k|15k|20k|30k|40k|50k|70k|100k|150k|200k|300k
            Preset Value: 70k
            Unit:Hz
        '''
        
        self.write(":SENSe1:BANDwidth:RESolution {}".format(freq))
        
    @property
    def center_frequency(self):
        '''This command gets the center frequency.
        
            Parameter: Value
            Description: Center Frequency
            Data Type: Numeric type (Real)
            Range: 100kHz to maximum frequency
            Preset Value: Depending on frequency option
            Unit: Hz
            Resolution: 500m
        '''
        
        freq = int(float(self.query(":SENSe1:FREQuency:CENTer?").split("\n")[0]))
        self._center_freq = freq
        return self._center_freq
    @center_frequency.setter
    def center_frequency(self,freq):
        '''This command sets the center frequency.
        
            Parameter: Value
            Description: Center Frequency
            Data Type: Numeric type (Real)
            Range: 100kHz to maximum frequency
            Preset Value: Depending on frequency option
            Unit: Hz
            Resolution: 500m
        '''
        
        self.write(":SENSe1:FREQuency:CENTer {}".format(freq))
    
    @property
    def start_frequency(self):
        '''This command sets the start frequency.
        
            Parameter: Value
            Description: Start frequency 
            Data Type: Numeric type (Real)
            Range: 100k to maximum frequency
            Preset Value: 100k
            Unit: Hz
            Resolution: 1
        '''
        
        freq = int(float(self.query(":SENSe1:FREQuency:start?").split("\n")[0]))
        self._start_freq = freq
        return self._start_freq
    @start_frequency.setter
    def start_frequency(self,freq):
        '''This command gets the start frequency.
        
            Parameter: Value
            Description: Start frequency 
            Data Type: Numeric type (Real)
            Range: 100k to maximum frequency
            Preset Value: 100k
            Unit: Hz
            Resolution: 1
        '''
        
        self.write(":SENSe1:FREQuency:start {}".format(freq))
    
    @property
    def stop_frequency(self):
        '''This command gets the stop frequency.
        
            Parameter: Value
            Description: Stop value
            Data Type: Numeric type (Real)
            Range: 100k ~ Maximum frequency
            Preset Value: Depending on frequency option
            Unit: Hz
            Resolution: 1
        '''
        
        freq = int(float(self.query(":SENSe1:FREQuency:stop?").split("\n")[0]))
        self._stop_freq = freq
        return self._stop_freq    
    @stop_frequency.setter
    def stop_frequency(self,freq):
        '''This command sets the stop frequency.

        Parameter: Value
        Description: Stop value
        Data Type: Numeric type (Real)
        Range: 100k ~ Maximum frequency
        Preset Value: Depending on frequency option
        Unit: Hz
        Resolution: 1
        '''
        self.write(":SENSe1:FREQuency:stop {}".format(freq))
    
    @property
    def span_frequency(self):
        '''This command gets the frequency span.
        
            Parameter:Value
            Description: Span frequency
            Data Type: Numeric type (Real)
            Range: 0 ~ Maximum Frequency -100k 
            Preset Value: Depending on frequency option
            Unit: Hz
            Resolution: 1
        '''
        freq = int(float(self.query(":SENSe1:FREQuency:span?").split("\n")[0]))
        self._span_freq = freq
        return self._span_freq        
    @span_frequency.setter
    def span_frequency(self,freq):
        '''This command sets the frequency span.
        
            Parameter:Value
            Description: Span frequency
            Data Type: Numeric type (Real)
            Range: 0 ~ Maximum Frequency -100k 
            Preset Value: Depending on frequency option
            Unit: Hz
            Resolution: 1
        '''
        self.write(":SENSe1:FREQuency:span {}".format(freq))
    
    @property
    def sweep_points(self):
        '''This command gets the number of measurement points.
        
            Parameter: Value
            Description: Number of measurement points
            Data Type: Numeric type (Integer)
            Range: 2 ~ 10001
            Preset Value: 201
            Resolution: 1
        '''
        
        points = int(self.query(":SENSe1:SWEep:POINts?").split("\n")[0])
        self._sweep_points = points
        return self._sweep_points
    @sweep_points.setter
    def sweep_points(self,points):
        '''This command sets the number of measurement points.
        
            Parameter: Value
            Description: Number of measurement points
            Data Type: Numeric type (Integer)
            Range: 2 ~ 10001
            Preset Value: 201
            Resolution: 1
        '''
        
        self.write(":SENSe1:SWEep:POINts {}".format(points))
        
    @property
    def average_points(self):
        '''This command gets the averaging factor.
        
            Parameter: Value
            Description: Averaging factor
            Data Type: Numeric type (Integer)
            Range: 1 ~ 999
            Preset Value: 16
            Unit: -
            Resolution: -
        '''
        
        avg = int(self.query(":SENSe1:AVERage:COUNt?").split("\n")[0])
        self._average_points = avg
        return self._average_points
    @average_points.setter
    def average_points(self,avg):
        '''This command sets the averaging factor.
        
            Parameter: Value
            Description: Averaging factor
            Data Type: Numeric type (Integer)
            Range: 1 ~ 999
            Preset Value: 16
            Unit: -
            Resolution: -
        '''
        self.write(":SENSe1:AVERage:COUNt {}".format(avg))