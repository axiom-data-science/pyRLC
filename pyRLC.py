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
from scipy.interpolate import interp1d


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

        # scandata(NscanX, Nsamp)
        # scandata - X dim - NscanX
        self.NscanX = sweep_header['scanlines_per_sweep_ext']
        # scandata - Y dim - Nsamp
        self.Nsamp = sweep_header['samples_per_scanline']
        self.scandata = np.zeros((self.NscanX,
                                  self.Nsamp))
        # scanline numbers used for reprojection - scan_i
        self.scan_i = np.zeros((self.NscanX,))

        scan = 0
        fpos1 = f.tell()
        while fpos1 < self.file_size:
            try:
                scan_header = self.read_scan_header(f)
            except:
                print('exception: {}, {}'.format(fpos1, self.file_size))

            if scan_header['utype'] == 0:
                seg_info = self.read_extended_segment_info(f)

                if (seg_info['number'] >= 0) & (seg_info['number'] <= self.NscanX):
                    if scan_header['compressed']:
                        fpos0 = f.tell()
                        self.scandata[scan, :] = self.read_compressed_scanline(f)
                        fpos1 = f.tell()
                    else:
                        self.scandata[scan, :] = np.frombuffer(f.read(self.NscanX), np.uint8)
                    self.scan_i[scan] = seg_info['number']
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

        Nscan = scan
        scandata = self.scandata[:Nscan, :]
        scan_i = self.scan_i[:Nscan]

        return (scandata, Nscan, scan_i) 

    def read_sweep_header(self, f):
        """Read and return individual sweep header, save Nscanx and Nsamp"""
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

    def read_compressed_scanline(self, f):
        """Decompress run-length compressed scanline"""
        scanline = np.zeros((self.Nsamp,))
        fpos0 = f.tell()
        samprem = self.file_size - fpos0
        samp = 0 

        while (samp < self.Nsamp) & (samp <= samprem):
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

    def interp_scandata(self):
        """Interpolate and return scan data"""
        interp_scandata = np.zeros((self.NscanX, self.Nsamp))
       
        # slow, but similar to MATLAB method
        # - slight differences because MATLAB backfills, this forward fills 
        for s in np.arange(0, self.Nsamp):
            f = interp1d(self.scan_i, self.scandata[:, s], 'nearest', fill_value='extrapolate')
            interp_scandata[:, s] = f(np.arange(self.NscanX))
        
        return interp_scandata


    def project_to_circular_image(self):
        """Project data to circular radar image."""
        imgsz = (2*self.Nsamp) - 1
        rdrimg = np.zeros((imgsz, imgsz)) - 1

        # X & Y image coordinates
        ix = np.tile(np.arange(imgsz), (imgsz, 1)) 
        iy = ix.T

        # X & Y coords
        rx = ix - self.Nsamp
        ry = iy - self.Nsamp

        # rotdeg to 0 because the radar is aligned to north
        rotdeg = 0

        coord_range = np.sqrt(rx*rx + ry*ry)
        az = np.arctan2(ry, rx) + rotdeg*np.pi/180.0
        az = az + 2*np.pi*(az <= 0)

        samp = np.round(coord_range) 
        scan = np.round(self.NscanX*az/(2*np.pi))
        scan = scan + self.NscanX*(scan <= 0) - 1

        for x in range(imgsz):
            for y in range(imgsz):
                # scan and sample coords
                sa = int(samp[y, x])

                # for given x, y pixel
                sc = int(scan[y, x]) 

                if sa < self.Nsamp:
                    rdrimg[y, x] = self.interp_scandata[sc, sa]

        return rdrimg

    
def read_uint8(f):
    return np.frombuffer(f.read(1), np.uint8)[0]

def read_uint32(f):
    return np.frombuffer(f.read(4), np.uint32)[0]

def read_uint64(f):
    return np.frombuffer(f.read(8), np.uint64)[0]
