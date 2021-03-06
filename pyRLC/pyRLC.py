#!python
# utf-8
"""Methods for reading .rlc files and converting to .png and .geotiff.

rlc File format (type 4):
- Header - 8 bytes
- Sweeps
    - Header - 7 bytes
    - Scan data - dependent on file
"""
import argparse
import os
import logging
import sys

import numpy as np
from matplotlib.image import imsave
from scipy.interpolate import interp1d

from PIL import Image


class RLCfile:
    """Methods for manipulating .rcl files."""
    NMEA_SIZE = 240

    def __init__(self, fpath):
        self.file_path = os.path.abspath(fpath)
        self.file_name = os.path.split(self.file_path)[-1]
        self.file_size = os.path.getsize(self.file_path)
        # file handle
        self.file = None
        # output file name
        self.output_file = None

        # Scan information read from header
        # rlc type
        self.rlctype = None
        # scandata(NscanX, Nsamp)
        self.scandata = None
        # X dimension in scandata
        self.NscanX = None
        # Y dimension in scandata
        self.Nsamp = None
        # Scanline numbers used in projection to circle for image creation
        self.scan_i = None
        # NMEA header
        self.NMEAhead = None

        # interpolated scandata to square for circle image
        self.interpolated_scandata = None

        self.logger = logging.getLogger()

    # context manager
    def __enter__(self):
        self.file = open(self.file_path, 'rb')
        return self

    def __exit__(self, exx_type, exc_val, exc_tb):
        if self.file:
            self.file.close()

    def read_header(self):
        """Read and return header data"""
        filesig = np.frombuffer(self.file.read(4), np.uint8)
        rlctype = np.frombuffer(self.file.read(1), np.uint8)
        three_zeros = np.frombuffer(self.file.read(3), np.uint8)

        self.logger.info('Read header')
        return filesig, rlctype, three_zeros

    def read_record_data_type4(self):
        """Read and return record data (RLC type 4)"""
        sweep_header = self.read_sweep_header()

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
        fpos1 = self.file.tell()
        while fpos1 < self.file_size:
            try:
                scan_header = self.read_scan_header()
            except:
                self.logger.exception('fpos1: {}, file_size: {}'.format(fpos1, self.file_size))

            if scan_header['utype'] == 0:
                seg_info = self.read_extended_segment_info()

                if (seg_info['number'] >= 0) & (seg_info['number'] <= self.NscanX):
                    if scan_header['compressed']:
                        fpos0 = self.file.tell()
                        self.scandata[scan, :] = self.read_compressed_scanline()
                        fpos1 = self.file.tell()
                        self.logger.debug('{} bytes read'.format(fpos1-fpos0))
                    else:
                        self.scandata[scan, :] = np.frombuffer(self.file.read(self.NscanX), np.uint8)
                    self.scan_i[scan] = seg_info['number']
                    scan += 1
                else:
                    # end of file
                    fpos1 = self.file_size + 1
            elif scan_header['utype'] > 0:
                # read NMEA header 
                self.NMEAhead = self.read_NMEA_header()
            else:
                # At end of file
                fpos1 = self.file_size + 1

        Nscan = scan
        scandata = self.scandata[:Nscan, :]
        scan_i = self.scan_i[:Nscan]

        return scandata, Nscan, scan_i

    def read_sweep_header(self):
        """Read and return individual sweep header, save Nscanx and Nsamp"""
        sweep_header = {
            'sample_rate': read_uint32(self.file),
            'samples_per_scanline': read_uint32(self.file),
            'scanlines_per_sweep': read_uint32(self.file),
            'scanlines_per_sweep_ext': read_uint32(self.file),
            'time_of_sweep_pc': read_uint32(self.file),
            'time_of_sweep_vp': read_uint64(self.file),
            'uspare': read_uint32(self.file),
        }
        return sweep_header

    def read_scan_header(self):
        """Read header for individual scan
        
        Notes:
        - bits 1-3 (sum) : utype (0: scan data; else: NMEA data)
        - bits 4-7 (sep) : urange
        - bit 8          : compressed (specifies if following is runlength encoded)
        
        """
        scanhead = read_uint8(self.file)

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

    def read_NMEA_header(self):
        """Read and return NMEA header"""
        return np.frombuffer(self.file.read(self.NMEA_SIZE), np.uint8, count=240)

    def read_extended_segment_info(self):
        """Read extended segment info"""
        number = read_uint32(self.file)
        time = read_uint32(self.file)

        if type(number) is not np.uint32:
            number = -1
            time = np.nan

        return {'number': number, 'time': time}

    def read_compressed_scanline(self):
        """Decompress run-length compressed scanline"""
        scanline = np.zeros((self.Nsamp,))
        fpos0 = self.file.tell()
        samprem = self.file_size - fpos0
        samp = 0 

        while (samp < self.Nsamp) & (samp <= samprem):
            backscat = read_uint8(self.file)
            if backscat == 255:
                repval = read_uint8(self.file) + 1
                backscat = read_uint8(self.file)

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
        radar_image = np.zeros((imgsz, imgsz)) - 1

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
        az = np.arctan2(ry, rx) + rotdeg*np.pi/180
        az = az + 2*np.pi*(az <= 0)

        # Add 1 here to compare with MATLAB version
        samp = (np.round(coord_range) + 1).astype(int)
        # Nscan, not NscanX
        scan = np.round(self.NscanX*az/(2.0*np.pi))
        scan = (scan + self.NscanX*(scan <= 0)).astype(int)

        for x in range(imgsz):
            for y in range(imgsz):
                # scan and sample coords
                # - for given x, y pixel
                # - change to 0-based index
                sa = samp[y, x] - 1
                sc = scan[y, x] - 1

                if sa < self.Nsamp:
                    radar_image[y, x] = self.interpolated_scandata[sc, sa]

        return radar_image

    def write_png(self,  output_file=None):
        """Write projected scandata to png"""
        if not output_file:
            # save to same dir as .rec file, just change to .png
            self.output_file = os.path.splitext(self.file_path)[0] + '.png'
        else:
            self.output_file = output_file
        imsave(self.output_file, self.radar_image, cmap='gray', vmin=0, vmax=255)

    def rec_to_png(self, output_file=None):
        """Reads .rec file and saves .png representation"""
        # Read file header and ensure it has correct signature
        fs, self.rlctype, sparebytes = self.read_header()
        if fs[0] != 206 or fs[1] != 4 or fs[2] != 90 or fs[3] != 27:
            raise Exception('file signature in header is incorrect: {}'.format(str(f) for f in fs))

        # Read data (assumes it is compressed type4)
        self.scandata, self.Nscan, self.scan_i = self.read_record_data_type4()

        # Interpolate to square
        self.interpolated_scandata = self.interp_scandata()

        # Reproject
        self.radar_image = self.project_to_circular_image()

        # Save image
        self.write_png(output_file)


class RadarOverlay:
    """Methods to add and manipulate overlays on radar images"""

    def __init__(self, radar_path, overlay_path, output_path):
        self.logger = logging.getLogger()

        self.radar_path = radar_path
        self.overlay_path = overlay_path
        self.output_path = output_path

    def add_overlay(self):
        """Adds overlay """
        with Image.open(self.radar_path) as src, Image.open(self.overlay_path) as overlay:
            self.logger.info('overlaying {} on {}'.format(self.overlay_path, self.radar_path))
            src.paste(overlay, (0, 0), overlay)
            src.save(self.output_path)
            self.logger.info('saved output combined image: {}'.format(self.output_path))


def parse_args(args):
    parser = argparse.ArgumentParser(description='Convert .rec file(s) to .png')
    parser.add_argument('in_file', help='comma seperated .rec files (e.g. f1.rec,f2.rec)')
    parser.add_argument('out_dir', help='output directory')
    parser.add_argument('-o', '--overlay', help='overlay image')

    return parser.parse_args()


def main():
    """Parse cli and iterate over input .rec files"""
    args = parse_args(sys.argv[1:])
    in_files = args.in_file.split(',')
    out_dir = args.out_dir
    overlay = args.overlay

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    for f in in_files:
        with RLCfile(f) as rlc_file:
            # Convert .rec file to .png
            filename = os.path.split(f)[-1]
            out_file = os.path.splitext(filename)[0] + '.png'
            output_file = os.path.join(out_dir, out_file)
            rlc_file.rec_to_png(output_file)

            if overlay:
                # Add overlay to .png
                overlay = RadarOverlay(rlc_file.output_file, overlay, rlc_file.output_file)
                overlay.add_overlay()


def read_uint8(f):
    return np.frombuffer(f.read(1), np.uint8)[0]


def read_uint32(f):
    return np.frombuffer(f.read(4), np.uint32)[0]


def read_uint64(f):
    return np.frombuffer(f.read(8), np.uint64)[0]


if __name__ == '__main__':
    main()
