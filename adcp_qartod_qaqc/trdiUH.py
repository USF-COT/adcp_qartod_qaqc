
from pycurrents.adcp.rdiraw import Multiread
import numpy as np
import numpy.ma as ma
from datetime import datetime

from adcp_qartod_qaqc.tests import (
    battery_flag_test,
    checksum_test,
    bit_test,
    orientation_test,
    sound_speed_test,
    correlation_magnitude_test,
    percent_good_test,
    current_speed_test,
    current_direction_test,
    horizontal_velocity_test,
    vertical_velocity_test,
    error_velocity_test,
    echo_intensity_test,
    range_drop_off_test,
    current_speed_gradient_test
)


class TRDIQAQC(object):
    """
    Performs QARTOD Data Quality Assurance and Control Tests
    on TRDI ADCP Data according to:
    "IOOS:
    Manual for Real-Time Quality Control of In-Situ Current Observations"

    Expects data in the same format as output by
    University of Hawaii's Multiread function

    By: Jeff Donovan <jdonovan@usf.edu> & Michael Lindemuth <mlindemu@usf.edu>
    University of South Florida
    College of Marine Science
    """

    @staticmethod
    def from_file(path, read_type, transducer_depth):
        """
        A convenience method to read in a file by path
        """
        m = Multiread(path, read_type)
        return TRDIQAQC(m.read(), transducer_depth)

    def __init__(self, multiread_data, transducer_depth):
        self.data = multiread_data
        self.transducer_depth = transducer_depth

        self.__read_timestamp()
        self.__read_configuration()
        self.__read_velocities()

        self.set_ensemble_bottom_stats()

    def __read_timestamp(self):
        yr = str(self.data.yearbase)
        mon = str(self.data.VL[0][2]).zfill(2)
        day = str(self.data.VL[0][3]).zfill(2)
        hr = str(self.data.VL[0][4]).zfill(2)
        min = str(self.data.VL[0][5]).zfill(2)
        sec = str(self.data.VL[0][6]).zfill(2)
        dt = '/'.join([yr, mon, day])
        tm = ':'.join([hr, min, sec])
        time_str = ' '.join([dt, tm])
        self.timestamp = datetime.strptime(time_str, '%Y/%m/%d %H:%M:%S')

    def __read_configuration(self):
        """
        Reads some shared configuration variables into
        more easily accessible class variables
        """

        self.NCells = self.data.NCells  # number of data bins
        self.raw_depth = self.data.dep  # uncorrected depth of each bin
        # actual depth of bin 1
        self.bin1depth = self.data.Bin1Dist + self.transducer_depth
        # vector containing the corrected depth of each bin
        self.depth = self.raw_depth + self.transducer_depth
        # TRDI ADCP operating frequency
        self.frequency = self.data.sysconfig['kHz']
        # beam angle is necessary for calculating last good bin
        self.beam_angle = self.data.sysconfig['angle']
        # orientation is either up of down and stored as a string
        if (self.data.sysconfig['up']):
            self.orientation = 'up'
        else:
            self.orientation = 'down'

        # this is in cm instead of m
        self.bin_size = self.data.FL.CellSize
        # this is in cm instead of m
        self.blank = self.data.FL.Blank
        # number of pings per ensemble (sample)
        self.pings = self.data.FL.NPings
        # lag is used to calculate distance to first good bin
        self.transmit_lag = self.data.FL.TransLag
        # pulse length is used to calculate distance to first good bin
        self.transmit_pulse = self.data.FL.Pulse
        self.TPMin = str(self.data.FL.TPP_min).zfill(2)  # ping minutes
        self.TPSec = str(self.data.FL.TPP_sec).zfill(2)  # ping seconds
        self.TP100 = str(self.data.FL.TPP_hun).zfill(2)  # ping hundreds
        # ping interval text string
        self.ping_interval = ':'.join([self.TPMin, self.TPSec, self.TP100])
        # magnetic declination applied to data, site specific, use  r input
        self.mag_declination = self.data.FL.EV / 10.0

    def __read_velocities(self):
        """
        extract velocity from the masked array
        Univ of Hawaii code converts the mm/s
        velocity data from the instrument to
        m/s. We work in cm/s so I will have to
        convert the velocities from m/s to cm/s
        """
        u_vel_masked = self.data.vel1
        u_vel_masked.mask = ma.nomask
        u_vel = u_vel_masked.compressed()
        self.u = u_vel * 100.

        v_vel_masked = self.data.vel2
        v_vel_masked.mask = ma.nomask
        v_vel = v_vel_masked.compressed()
        self.v = v_vel * 100.

        w_vel_masked = self.data.vel3
        w_vel_masked.mask = ma.nomask
        w_vel = w_vel_masked.compressed()
        self.w = w_vel * 100.

        err_vel_masked = self.data.vel4
        err_vel_masked.mask = ma.nomask
        err_vel = err_vel_masked.compressed()
        self.ev = err_vel * 100.
        self.ev = self.ev

        self.z = self.u + 1j + self.v
        self.current_speed = abs(self.z)
        self.current_direction = np.arctan2(self.z.real, self.z.imag)*180/np.pi
        self.current_direction = self.current_direction

    def set_ensemble_bottom_stats(self, tolerance=30):
        self.ensemble_bottom_stats = []
        for ensemble in self.data.amp:
            bottom_bin = 1
            for bin_prev, bin_curr in zip(ensemble, ensemble[1:]):
                bin_flag_count = 0
                for beam_prev, beam_curr in zip(bin_prev, bin_curr):
                    bin_diff = int(beam_curr) - int(beam_prev)
                    if bin_diff > tolerance:
                        bin_flag_count += 1

                if bin_flag_count < 2:
                    bottom_bin += 1
                else:
                    break

            bottom_stats = {}
            bottom_stats['bottom_bin'] = bottom_bin
            bottom_stats['range_to_bottom'] = (
                (bottom_stats['bottom_bin'] * (self.NCells / 100)) +
                self.bin1depth)
            bottom_stats['side_lobe_start'] = (
                (int(np.cos(self.data.sysconfig['angle'] * (np.pi/180.)) *
                 bottom_stats['range_to_bottom'])))
            bottom_stats['last_good_bin'] = bottom_stats['side_lobe_start'] - 1
            bottom_stats['last_good_counter'] = (
                bottom_stats['last_good_bin'] - 1)
            self.ensemble_bottom_stats.append(bottom_stats)

    def battery_flag(self):
        """
        QARTOD Test #1 Strongly Recommended
        battery flag test can not be performed
        if the GOES header f_code=G and the data can be
        decoded, then sample is good
        """
        return battery_flag_test(self.data)

    def checksum_flag(self):
        """
        QARTOD Test #2 Required
        Checksum test can not be performed because
        data is converted to/from Psuedo-ASCII during the
        GOES transmission. If the GOES header f_code=G
        and the data can be decoded, then sample is good
        """
        return checksum_test(self.data)

    def bit_flag(self):
        """
        Not an official QARTOD test.  Checks special TRDI bit flag.
        """
        bit_flags = []
        for ensemble in self.data.VL:
            bit_flags.append(bit_test(ensemble))

        return bit_flags

    def orientation_flags(self, max_pitch=20, max_roll=20):
        """
        QARTOD Test #3 Required: orientation (pitch and roll) tests
        """
        return orientation_test(self.data.pitch, self.data.roll,
                                max_pitch, max_roll)

    def sound_speed_flags(self, sound_speed_min=1400, sound_speed_max=1600):
        """
        QARTOD Test 4 Required: Sound speed test
        """
        ssval_flags = []
        for ensemble_ssv in self.data.VL:
            ssval_flags.append(
                sound_speed_test(ensemble_ssv,
                                 sound_speed_min, sound_speed_max))

        return ssval_flags

    # NOTE: QARTOD Tests 5, 6, and 7 cannot be performed on TRDI ADCP

    def correlation_magnitude_flags(self,
                                    good_tolerance=115,
                                    questionable_tolerance=64):
        """
        QARTOD Test #8 Strongly Recommended
        correlation magnitude test
        """
        ensemble_flags = []
        for i in xrange(0, len(self.data.cor1)):
            correlation = zip(self.data.cor1[i], self.data.cor2[i],
                              self.data.cor3[i], self.data.cor4[i])
            correlation_flags = correlation_magnitude_test(
                correlation,
                good_tolerance,
                questionable_tolerance
            )
            ensemble_flags.append(correlation_flags)

        return ensemble_flags

    def percent_good_flags(self, percent_good=21, percent_bad=17):
        """
        QARTOD Test #9 Required
        percent good test
        default limits on this test are > 21%, good: < 17%, bad
        limits derived from TRDI Spreadsheed based on our instruments and setup
        """
        ensemble_flags = []
        for pg3, pg4, bottom_stats in zip(self.data.pg3, self.data.pg4,
                                          self.ensemble_bottom_stats):
            ensemble_flags.append(
                percent_good_test(pg3, pg4, bottom_stats['last_good_counter'],
                                  percent_good, percent_bad)
            )

        return ensemble_flags

    def current_speed_flags(self, max_speed=150):
        """
        QARTOD Test #10 Required
        current speed test
        150 cm/s is the West Florida Shelf limit.  Adjust as necessary
        """

        last_good_counter = self.ensemble_bottom_stats[0]['last_good_counter']
        return current_speed_test(self.current_speed,
                                  last_good_counter)

    def current_direction_flags(self):
        """
        QARTOD Test #11 Required
        current direction test
        Verifies that the current direction is less than 360 degrees
        Negative values are made positive by adding 360 to them
        """

        last_good_bin = self.ensemble_bottom_stats[0]['last_good_counter']
        return current_direction_test(self.current_direction,
                                      last_good_bin)

    def horizontal_velocity_flags(self,
                                  max_u_vel=150, max_v_vel=150):
        """
        QARTOD Test #12 Required
        horizontal velocity test
        If EITHER u or v greater than 150 cm/s, velocity is bad
        150 cm/s is a local WFS limit
        """

        last_good_bin = self.ensemble_bottom_stats[0]['last_good_counter']
        return horizontal_velocity_test(self.u, self.v,
                                        last_good_bin,
                                        max_u_vel, max_v_vel)

    def vertical_velocity_flags(self, max_w_velocity=15):
        """
        QARTOD Test #13 Strongly Recommended
        vertical velocity test
        if w greater than 15 cm/s (10% MAX speed from TRDI), w is bad
        """

        last_good_bin = self.ensemble_bottom_stats[0]['last_good_counter']
        return vertical_velocity_test(self.w, last_good_bin,
                                      max_w_velocity)

    def error_velocity_flags(self,
                             questionable_error_velocity=2.6,
                             bad_error_velocity=5.2):
        """
        QARTOD Test #14 Required
        error velocity test
        err < 2.6 cm/s, good: err > 5.2 cm/s, bad
        limits derived from TRDI Spreadsheet
        based on our instruments and setup
        """

        last_good_bin = self.ensemble_bottom_stats[0]['last_good_counter']
        return error_velocity_test(self.ev, last_good_bin,
                                   questionable_error_velocity,
                                   bad_error_velocity)

    def echo_intensity_flags(self, tolerance=2):
        """
        QARTOD Test #15 Required
        echo intensity test
        """

        ensemble_flags = []
        for ensemble in self.data.amp:
            ensemble_flags.append(
                echo_intensity_test(ensemble, tolerance)
            )

        return ensemble_flags

    def range_drop_off_flags(self, drop_off_limit=60):
        """
        QARTOD Test #16 Strongly Recommended
        range drop-off test
        Range Limit set to 60 as recommended in QARTOD spreadsheet (CO-OPS).
        The QARTOD recommended cut-off is 30. ???
        """

        ensemble_flags = []
        for ensemble in self.data.amp:
            ensemble_flags.append(
                range_drop_off_test(ensemble, drop_off_limit)
            )

        return ensemble_flags

    def current_speed_gradient_flags(self, tolerance=6):
        """
        QARTOD Test #17 Strongly Recommended
        current speed gradient test
        """

        last_good_bin = self.ensemble_bottom_stats[0]['last_good_counter']
        return current_speed_gradient_test(self.current_speed,
                                           last_good_bin,
                                           tolerance)

import sys
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_path", help="Path of PD0 file to parse")
    parser.add_argument("read_type", help="""
    Multiread read type (e.g., wh).  See Mutiread documentation for details
    """)
    parser.add_argument("transducer_height",
                        type=float, help="Depth of ADCP transducer")
    args = parser.parse_args()

    trdi_qaqc = TRDIQAQC.from_file(args.input_path,
                                   args.read_type, args.transducer_height)
    print 'Ensemble Bottom Bins: %s' % (trdi_qaqc.ensemble_bottom_stats,)

    return 1

if __name__ == '__main__':
    sys.exit(main())
