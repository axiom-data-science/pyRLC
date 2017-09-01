import pyRLC
import numpy as np


TEST_FILE = 'tests/test.rec'


def test_read_header():
    rlc_file = pyRLC.RLCfile(TEST_FILE)
    
    with open(TEST_FILE, 'rb') as f:
        (filesig, rlctype, sparebytes) = rlc_file.read_header(f)

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
    rlc_file = pyRLC.RLCfile(TEST_FILE)

    with open(TEST_FILE, 'rb') as f:
        # skip header 
        f.seek(8)

        sweep_header = rlc_file.read_sweep_header(f)
        assert sweep_header['sample_rate'] == 54000000 
        assert sweep_header['samples_per_scanline'] == 512 
        assert sweep_header['scanlines_per_sweep'] == 1024
        assert sweep_header['scanlines_per_sweep_ext'] == 4096 
        assert sweep_header['time_of_sweep_pc'] == 1973129824 
        assert sweep_header['time_of_sweep_vp'] == 18313917820654764054 
        assert sweep_header['uspare'] == 0


def test_read_scan_header():
    rlc_file = pyRLC.RLCfile(TEST_FILE)

    with open(TEST_FILE, 'rb') as f:
        # skip header and sweep header
        f.seek(40)

        scan_header = rlc_file.read_scan_header(f)
        assert scan_header['utype'] == 1
        assert np.array_equal(scan_header['urange'], np.array([0, 0, 0, 0]))
        assert scan_header['compressed'] == 0 


def test_read_NMEA_header():
    rlc_file = pyRLC.RLCfile(TEST_FILE)

    with open(TEST_FILE, 'rb') as f:
        # skip header, sweep header, and scan head
        f.seek(41)

        NMEA_header = rlc_file.read_NMEA_header(f)
        assert NMEA_header[0] == 35
        assert NMEA_header[1] == 0
        assert NMEA_header[8] == 180
        assert NMEA_header[9] == 14
        assert NMEA_header[-1] == 0
        assert NMEA_header.sum() == 906


def test_read_scanhead2():
    rlc_file = pyRLC.RLCfile(TEST_FILE)

    with open(TEST_FILE, 'rb') as f:
        # skip header, sweep header, scan head, NMEA header
        f.seek(281)

        scanhead = rlc_file.read_scan_header(f)
        assert scanhead['utype'] == 0
        assert np.array_equal(scanhead['urange'], np.array([1, 0, 1, 0]))
        assert scanhead['compressed'] == 1 


def test_read_extended_segment_info():
    rlc_file = pyRLC.RLCfile(TEST_FILE)

    with open(TEST_FILE, 'rb') as f:
        # skip
        f.seek(282)

        extended_info = rlc_file.read_extended_segment_info(f)
        assert extended_info['number'] == 0
        assert extended_info['time'] == 4264016289 


def test_read_compressed_scanline():
    rlc_file = pyRLC.RLCfile(TEST_FILE)

    with open(TEST_FILE, 'rb') as f:
        # Read header and save required info
        f.seek(8) 

        sweep_header = rlc_file.read_sweep_header(f)
        # scandata(NscanX, Nsamp)
        # scandata - X dim - NscanX
        rlc_file.NscanX = sweep_header['scanlines_per_sweep_ext']
        # scandata - Y dim - Nsamp
        rlc_file.Nsamp = sweep_header['samples_per_scanline']
        rlc_file.scandata = np.zeros((rlc_file.NscanX,
                                      rlc_file.Nsamp))
        # scanline numbers used for reprojection - scan_i
        rlc_file.scanline_nums = np.zeros((rlc_file.NscanX,))

        # skip to my lou
        f.seek(290)
        scandata = rlc_file.read_compressed_scanline(f)
        assert scandata[0] == 5
        assert scandata[-1] == 0
        assert scandata.sum() == 6522


def test_read_compressed_scanline2():
    rlc_file = pyRLC.RLCfile(TEST_FILE)

    with open(TEST_FILE, 'rb') as f:
        # Read header and save required info
        f.seek(8) 
        sweep_header = rlc_file.read_sweep_header(f)
        # scandata(NscanX, Nsamp)
        # scandata - X dim - NscanX
        rlc_file.NscanX = sweep_header['scanlines_per_sweep_ext']
        # scandata - Y dim - Nsamp
        rlc_file.Nsamp = sweep_header['samples_per_scanline']
        rlc_file.scandata = np.zeros((rlc_file.NscanX,
                                      rlc_file.Nsamp))
        # scanline numbers used for reprojection - scan_i
        rlc_file.scanline_nums = np.zeros((rlc_file.NscanX,))

        # skip to my lou
        f.seek(799)

        scandata = rlc_file.read_compressed_scanline(f)
        assert scandata[0] == 168 
        assert scandata[-1] == 10 
        assert scandata.sum() == 7213 


def test_read_rec_data():
    rlc_file = pyRLC.RLCfile(TEST_FILE)

    with open(TEST_FILE, 'rb') as f:
        # Read header and save required info
        f.seek(8) 
        (scandata, Nscan, scan_i) = rlc_file.read_record_data_type4(f)
        assert Nscan == 2453
        assert scan_i.sum() == 5043250
        assert scandata.sum().sum() == 16223096


def test_interp_scandata():
    rlc_file = pyRLC.RLCfile(TEST_FILE)

    with open(TEST_FILE, 'rb') as f:
        # Read header and save required info
        f.seek(8) 
        (rlc_file.scandata, rlc_file.Nscan, rlc_file.scan_i) = rlc_file.read_record_data_type4(f)

        interp_scandata = rlc_file.interp_scandata(f) 
        assert interp_scandata.sum().sum() == 27092281 
