import pyRLC


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
