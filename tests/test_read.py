# Numbers for assert values derived from MATLAB implementation
import os
import numpy as np

from scipy.io import loadmat
from scipy.misc import imread

import pyRLC

RADAR_FILE = 'data/test.rec'
MATLAB_FILE = 'data/interp_scandata.mat'
TEST_IMAGE = 'data/test_image.png'
REF_IMAGE = 'data/reference_image.png'


def test_read_header():
    with pyRLC.RLCfile(RADAR_FILE) as f:
        (filesig, rlctype, sparebytes) = f.read_header()

        # 4-byte signature (206 4 90 27)
        assert filesig[0] == 206
        assert filesig[1] == 4
        assert filesig[2] == 90
        assert filesig[3] == 27

        # RLC Type
        assert rlctype == 6

        # (0 0 0) ¯\_(ツ)_/¯
        assert sparebytes[0] == 0
        assert sparebytes[1] == 0
        assert sparebytes[2] == 0


def test_read_sweep_header():
    with pyRLC.RLCfile(RADAR_FILE) as f:
        # skip header 
        f.file.seek(8)

        sweep_header = f.read_sweep_header()
        assert sweep_header['sample_rate'] == 54000000 
        assert sweep_header['samples_per_scanline'] == 512 
        assert sweep_header['scanlines_per_sweep'] == 1024
        assert sweep_header['scanlines_per_sweep_ext'] == 4096 
        assert sweep_header['time_of_sweep_pc'] == 1973129824 
        assert sweep_header['time_of_sweep_vp'] == 18313917820654764054 
        assert sweep_header['uspare'] == 0

def test_read_scan_header():
    with pyRLC.RLCfile(RADAR_FILE) as f:
        # skip header and sweep header
        f.file.seek(40)

        scan_header = f.read_scan_header()
        assert scan_header['utype'] == 1
        assert np.array_equal(scan_header['urange'], np.array([0, 0, 0, 0]))
        assert scan_header['compressed'] == 0 


def test_read_NMEA_header():
    with pyRLC.RLCfile(RADAR_FILE) as f:
        # skip header, sweep header, and scan head
        f.file.seek(41)

        NMEA_header = f.read_NMEA_header()
        assert NMEA_header[0] == 35
        assert NMEA_header[1] == 0
        assert NMEA_header[8] == 180
        assert NMEA_header[9] == 14
        assert NMEA_header[-1] == 0
        assert NMEA_header.sum() == 906


def test_read_scanhead2():
    with pyRLC.RLCfile(RADAR_FILE) as f:
        # skip header, sweep header, scan head, NMEA header
        f.file.seek(281)

        scanhead = f.read_scan_header()
        assert scanhead['utype'] == 0
        assert np.array_equal(scanhead['urange'], np.array([1, 0, 1, 0]))
        assert scanhead['compressed'] == 1 


def test_read_extended_segment_info():
    with pyRLC.RLCfile(RADAR_FILE) as f:
        # skip
        f.file.seek(282)

        extended_info = f.read_extended_segment_info()
        assert extended_info['number'] == 0
        assert extended_info['time'] == 4264016289 


def test_read_compressed_scanline():
    with pyRLC.RLCfile(RADAR_FILE) as f:
        # Read header and save required info
        f.file.seek(8)

        sweep_header = f.read_sweep_header()
        # scandata(NscanX, Nsamp)
        # scandata - X dim - NscanX
        f.NscanX = sweep_header['scanlines_per_sweep_ext']
        # scandata - Y dim - Nsamp
        f.Nsamp = sweep_header['samples_per_scanline']
        f.scandata = np.zeros((f.NscanX, f.Nsamp))
        # scanline numbers used for reprojection - scan_i
        f.scanline_nums = np.zeros((f.NscanX,))

        # skip to correct section
        f.file.seek(290)
        scandata = f.read_compressed_scanline()
        assert scandata[0] == 5
        assert scandata[-1] == 0
        assert scandata.sum() == 6522


def test_read_compressed_scanline2():
    with pyRLC.RLCfile(RADAR_FILE) as f:
        # Read header and save required info
        f.file.seek(8)
        sweep_header = f.read_sweep_header()
        # scandata(NscanX, Nsamp)
        # scandata - X dim - NscanX
        f.NscanX = sweep_header['scanlines_per_sweep_ext']
        # scandata - Y dim - Nsamp
        f.Nsamp = sweep_header['samples_per_scanline']
        f.scandata = np.zeros((f.NscanX, f.Nsamp))
        # scanline numbers used for reprojection - scan_i
        f.scanline_nums = np.zeros((f.NscanX,))

        # skip to correct section
        f.file.seek(799)

        scandata = f.read_compressed_scanline()
        assert scandata[0] == 168 
        assert scandata[-1] == 10 
        assert scandata.sum() == 7213 


def test_read_record_data_type4():
    with pyRLC.RLCfile(RADAR_FILE) as f:
        # Read header and save required nfo
        f.file.seek(8)
        (scandata, Nscan, scan_i) = f.read_record_data_type4()
        assert Nscan == 2453
        assert scan_i.sum() == 5043250
        assert scandata.sum().sum() == 16223096


def test_interp_scandata():
    with pyRLC.RLCfile(RADAR_FILE) as f:
        # Read header and save required info
        f.file.seek(8)
        (f.scandata, f.Nscan, f.scan_i) = f.read_record_data_type4()

    interp_scandata = f.interp_scandata()
    # NOTE:  MATLAB nearest neighbor algorithm fills differently than Numpy (Scipy) implementations
    # O P M
    # 1 1 1
    # - 1 2
    # 2 2 2
    # - 2 3
    # - 2 3
    # 3 3 3
    # The assertion here is for the Scipy version, in MATLAB the value is 27088914
    assert interp_scandata.sum() == 27092281


def test_reprojection():
    with pyRLC.RLCfile(RADAR_FILE) as f:
        # Read header and save required info
        f.file.seek(8)
        (f.scandata, f.Nscan, f.scan_i) = f.read_record_data_type4()

    f.interpolated_scandata = f.interp_scandata()
    radar_image = f.project_to_circular_image()

    # MATLAB value 10583729
    # Python value 10583770
    assert radar_image.sum() == 10583770


def test_reprojection_matlab():
    """Test reproject with MATLAB interpolated values"""
    with pyRLC.RLCfile(RADAR_FILE) as f:
        # Read header and save required info
        f.file.seek(8)
        f.scandata, f.Nscan, f.scan_i = f.read_record_data_type4()

    # Read scan data array from MATLAB for test
    matlab_scandata = loadmat(MATLAB_FILE)['scandata']
    f.interpolated_scandata = matlab_scandata
    radar_image = f.project_to_circular_image()

    # MATLAB value 10583729
    # Python value 10583770
    assert radar_image.sum() == 10583729


def test_radar_image():
    with pyRLC.RLCfile(RADAR_FILE) as f:
        # Read header and save required info
        f.file.seek(8)
        f.scandata, f.Nscan, f.scan_i = f.read_record_data_type4()

    # Read scan data array from MATLAB for test to create same image
    matlab_scandata = loadmat(MATLAB_FILE)['scandata']
    f.interpolated_scandata = matlab_scandata
    f.radar_image = f.project_to_circular_image()

    # Save image
    f.write_png(TEST_IMAGE)

    # Compare image with reference
    ref_image = imread(REF_IMAGE)
    test_image = imread(TEST_IMAGE)
    assert np.array_equal(ref_image, test_image)
