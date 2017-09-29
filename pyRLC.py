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
import logging

from matplotlib.image import imsave
import numpy as np
from scipy.interpolate import interp1d


class RLCfile:
    """Methods for manipulating .rcl files."""
    NMEA_SIZE = 240

    def __init__(self, fpath):
        self.file_path = os.path.abspath(fpath) 
        self.file_size = os.path.getsize(self.file_path)

        # Scan information read from header
        # rlc type
        self.rlctype = None
        # scandata(NscanX, Nsamp)
        self.scandata = None
        # X dimension in scandata
        self.NscanX = None
        # Y dimension in scandata
        self.Nsamp = None
        # Scanline numbers used in reprojection
        self.scan_i = None
        # NMEA header
        self.NMEAhead = None

#        with open(file_path, 'rb') as f:
#            (filesig, rlctype, _) = read_header(f)
        # logger
        self.logger = logging.getLogger()

    def read_header(self, f):
        """Read and return header data"""
        filesig = np.frombuffer(f.read(4), np.uint8)
        rlctype = np.frombuffer(f.read(1), np.uint8)
        three_zeros = np.frombuffer(f.read(3), np.uint8)

        self.logger.info('Read header')
        return filesig, rlctype, three_zeros

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

        self.logger.debug('NscanX: {}, Nsamp: {}'.format(self.NscanX, self.Nsamp))

        scan = 0
        fpos1 = f.tell()
        while fpos1 < self.file_size:
            try:
                scan_header = self.read_scan_header(f)
            except:
                self.logger.exception('fpos1: {}, file_size: {}'.format(fpos1, self.file_size))

            if scan_header['utype'] == 0:
                seg_info = self.read_extended_segment_info(f)

                if (seg_info['number'] >= 0) & (seg_info['number'] <= self.NscanX):
                    if scan_header['compressed']:
                        fpos0 = f.tell()
                        self.scandata[scan, :] = self.read_compressed_scanline(f)
                        fpos1 = f.tell()
                        self.logger.debug('{} bytes read'.format(fpos1-fpos0))
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

        return scandata, Nscan, scan_i

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
# Use?
#    def read_scanline(self, f):
#        """Read scanline"""
#        if self.scanhead['compressed'] == 0:
#            self.scandata[scan, :] = to_int(f.read(Nsamp))
# not used
#            self.scanlines[scan] = 1
#        else:
#            fpos0 = f.tell()
#            self.scandata[scan, :] = self.decompress_scanline(f)
# not used
#            self.scanlines[scan] = 1
#            fpos1 = f.tell()

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
        """Project data to circular projection"""
        imgsz = (2*self.Nsamp) - 1
        rdrimg = np.zeros((imgsz, imgsz)) - 1

        # X & Y image coordinates
        ix = np.tile(np.arange(1, imgsz+1), (imgsz, 1))
        iy = ix.T

        # X & Y coords
        rx = ix - self.Nsamp
        ry = iy - self.Nsamp

        # rotdeg to 0 because the radar is aligned to north
        rotdeg = 0

        # range and azimuth coordinates
        coord_range = np.sqrt(rx*rx + ry*ry)
        # relative to radar, which is 0 because radar is aligned to North
        #### - atan is slightly different :/
        az = np.arctan2(ry, rx) + rotdeg*np.pi/180
        az = az + 2*np.pi*(az <= 0)

        samp = np.round(coord_range) + 1
        # Nscan, not NscanX
        scan = np.round(self.NscanX*az/(2.0*np.pi))
        scan = scan + self.NscanX*(scan <= 0)

        for x in range(imgsz):
            for y in range(imgsz):
                # scan and sample coords
                # - for given x, y pixel
                # - change to 0-based index
                sa = int(samp[y, x]) - 1
                sc = int(scan[y, x]) - 1

                if sa < self.Nsamp:
                    rdrimg[y, x] = self.interp_scandata[sc, sa]


        return rdrimg

    def write_png(self, output_file):
        """Write projected scandata to png"""
        # Cast to uint8 to match MATLAB
#        rdrimg = self.rdrimg.astype('u1')
#        rdrimg[rdrimg == 255] = 0
        rdrimg = self.rdrimg
        imsave(output_file, rdrimg, cmap='gray')

    def do_this(self):
        """Reads .rec file and saves .png representation"""

def main():
    """Parse cli and iterate over input .rec files"""
    import os
    import argparse

    parser = argparse.ArgumentParser(description='Convert .rec file(s) to .png')
    parser.add_argument('in_file', help='comma seperated .rec files (e.g. f1.rec,f2.rec)')
    parser.add_argument('out_dir', help='output directory')
    args = parser.parse_args()

    in_files = args.in_file.split(',')
    out_dir = args.out_dir

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    for f in in_files:
        pass


def read_uint8(f):
    return np.frombuffer(f.read(1), np.uint8)[0]


def read_uint32(f):
    return np.frombuffer(f.read(4), np.uint32)[0]


def read_uint64(f):
    return np.frombuffer(f.read(8), np.uint64)[0]


if __name__ == '__main__':
    main()
