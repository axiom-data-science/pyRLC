#!python
# utf-8
"""Methods for reading .rlc files and converting to .png and .geotiff.

rlc File format (type 4):
- Header - 8 bytes
- Sweeps
    - Header - 7 bytes
    - Scan data - dependent on file
"""
import os

import numpy as np


class RLCfile:
    """All methods for manipulating .rcl files."""
    NMEA_SIZE = 240

    def __init__(self, fpath):
        self.file_path = os.path.abspath(fpath) 
        self.file_size = os.path.getsize(self.file_path) 

#        with open(file_path, 'rb') as f:
#            (filesig, rlctype, _) = read_header(f)

    def read_header(self, f):
        """Read and return header data"""
        filesig = np.frombuffer(f.read(4), np.uint8)
        rlctype = np.frombuffer(f.read(1), np.uint8)
        three_zeros = np.frombuffer(f.read(3), np.uint8)

        return (filesig, rlctype, three_zeros)


    def read_record_data_type4(self, f):
        """Read and return record data (RLC type 4)"""
        sweep_header = self.read_sweep_header(f)
        # scandata - X dim - NscanX
        self.samples_per_scan = sweep_header['scanlines_per_sweep_ext']
        # scandata - Y dim - Nsamp
        self.samples_per_sweep = sweep_header['samples_per_scanline']
        self.scandata = np.zeros((self.samples_per_scan,
                                  self.samples_per_sweep))
        # scanline numbers used for reprojection - scan_i
        self.scanline_nums = np.zeros((self.samples_per_scan,))

        scan = 0
        fpos1 = f.tell()
        while fpos1 <= self.file_size:
            scan_header = self.read_scan_header(f)

            if scan_header['utype'] == 0:
                seg_info = self.read_extended_segment_info(f)

                if (seg_info['number'] >= 0) & (seg_info['number'] <= self.samples_per_sweep):
                    if scan_header['compressed']:
                        self.scandata[scan, :] = self.read_compressed_scanline(f)
                    else:
                        fpos0 = f.tell()
                        self.scandata[scan, :] = np.frombuffer(f.read(self.samples_per_sweep), np.uint8)
                        fpos1 = f.tell()
                    self.scanline_nums[scan] = seg_info['number']
                    scan += 1
                else:
                    # end of file
                    fpos1 = self.file_size + 1
            elif scan_header['utype'] > 0:
                # read NMEA header 
                self.NMEAhead = self.read_NMEA_header(f)
            else:
                # At end of file
                fpos1 = self.file_size + 1

        Nscan = scan - 1
        scandata = self.scandata[:Nscan, :]
        scan_i = self.scanline_nums[:Nscan]

        return (scandata, Nscan, scan_i) 

    def read_sweep_header(self, f):
        """Read and return individual sweep header"""
        sweep_header = {}
        sweep_header['sample_rate'] = read_uint32(f)
        sweep_header['samples_per_scanline'] = read_uint32(f) 
        sweep_header['scanlines_per_sweep'] = read_uint32(f) 
        sweep_header['scanlines_per_sweep_ext'] = read_uint32(f) 
        sweep_header['time_of_sweep_pc'] = read_uint32(f) 
        sweep_header['time_of_sweep_vp'] = read_uint64(f) 
        sweep_header['uspare'] = read_uint32(f) 

        return sweep_header

    def read_scan_header(self, f):
        """Read header for individual scan
        
        Notes:
        - bits 1-3 (sum) : utype (0: scan data; else: NMEA data)
        - bits 4-7 (sep) : urange
        - bit 8          : compressed (specifies if following is runlength encoded)
        
        """
        scanhead = read_uint8(f) 

        if scanhead > 0:
            rsi = np.zeros((8,))
            b = 0
            while scanhead > 0:
                rsi[b] = scanhead % 2
                scanhead = scanhead // 2
                b += 1

            utype = rsi[:3].sum()
            urange = rsi[3:7]
            compressed = rsi[-1]
        else:
            utype = -1
            urange = -1
            compressed = -1
        
        return {'utype': utype, 'urange': urange, 'compressed': compressed}

    def read_NMEA_header(self, f):
        """Read and return NMEA header"""
        return np.frombuffer(f.read(self.NMEA_SIZE), np.uint8, count=240)

    def read_extended_segment_info(self, f):
        """Read extended segment info"""
        number = read_uint32(f) 
        time = read_uint32(f) 

        if type(number) is not np.uint32:
            number = -1
            time = np.nan

        return {'number': number, 'time': time}

    def read_scanline(self, f):
        """Read scanline"""
        if self.scanhead['compressed'] == 0:
            self.scandata[scan, :] = to_int(f.read(Nsamp)) 
#            self.scanlines[scan] = 1
        else:
            fpos0 = f.tell()
            self.scandata[scan, :] = self.decompress_scanline(f)
#            self.scanlines[scan] = 1
            fpos1 = f.tell()

    def _read_compressed_scanline(self, f):
        """Read run-length encoded scanline 
        
        Notes:
        - Corresponds to Xenex/RTI's "RLC_New"
        """
        self.fpos1 = f.tell()
        self.scandata[scan, :] = self.decompress_scanline(f)
        self.scanlines[scan] = 1
        self.fpos2 = f.tell()

    def read_compressed_scanline(self, f):
        """Decompress run-length compressed scanline"""
        Nsamp = self.samples_per_sweep
        scanline = np.zeros((self.samples_per_sweep,))
        fpos0 = f.tell()
        samprem = self.file_size - fpos0
        samp = 0 

        while (samp < Nsamp) & (samp <= samprem):
            backscat = read_uint8(f)
            if backscat == 255:
                repval = read_uint8(f) + 1
                backscat = read_uint8(f)

                scanline[samp:samp+repval] = backscat
                samp += repval
            else:
                scanline[samp] = backscat
                samp += 1

        return scanline


def read_uint8(f):
    return np.frombuffer(f.read(1), np.uint8)[0]

def read_uint32(f):
    return np.frombuffer(f.read(4), np.uint32)[0]

def read_uint64(f):
    return np.frombuffer(f.read(8), np.uint64)[0]
