import unittest

from trdi_adcp_readers.readers import (
    read_PD15_file
)

from adcp_qartod_qaqc.trdi import (
    TRDIQAQC
)

from adcp_qartod_qaqc.tests import (
    ADCP_FLAGS
)


class TestBaseTests(unittest.TestCase):
    pass


class TestTRDIQAQC(unittest.TestCase):

    def setUp(self):
        self.trdi_data = read_PD15_file('./test_data/1407B0B6', header_lines=2)
        self.qaqc = TRDIQAQC(self.trdi_data, transducer_depth=104)
        print self.trdi_data

    def test_bottom_stats(self):
        self.assertEqual(20, self.qaqc.bottom_stats['last_good_bin'])

    def test_battery_flag(self):
        self.assertEqual(ADCP_FLAGS['no_test'], self.qaqc.battery_flag())

    def test_checksum(self):
        self.assertEqual(ADCP_FLAGS['good'], self.qaqc.checksum_flag())

if __name__ == '__main__':
    unittest.main()
