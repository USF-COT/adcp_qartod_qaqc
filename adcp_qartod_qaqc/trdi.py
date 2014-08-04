
import math
from itertools import imap

from operator import itemgetter

from adcp_qartod_qaqc.tests import (
    battery_flag_test,
    checksum_test,
    bit_test,
    orientation_test,
    sound_speed_test,
    noise_floor_test,
    signal_strength_test,
    signal_to_noise_test,
    correlation_magnitude_test,
    percent_good_test,
    current_speed_test,
    current_direction_test,
    horizontal_velocity_test,
    vertical_velocity_test,
    error_velocity_test,
    stuck_sensor_test,
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

    Works on 1 ensemble at a time

    By: Jeff Donovan <jdonovan@usf.edu> & Michael Lindemuth <mlindemu@usf.edu>
    University of South Florida
    College of Marine Science
    """

    def __init__(self, data, transducer_depth=None):
        self.data = data

        if transducer_depth is not None:
            self.transducer_depth = transducer_depth
        else:
            self.transducer_depth = (
                self.data['variable_leader']['depth_of_transducer']
            )

        self.__read_velocities()
        self.__calc_bottom_stats()

    def __read_velocities(self):
        """
        Calculates z magnitude and direction given
        east (u), north(v), and vertical(w) velocities
        """

        self.u = []
        self.v = []
        self.w = []
        self.z = []
        self.current_speed = []
        self.current_direction = []
        for bin in self.data['velocity']['data']:
            u = bin[0]
            v = bin[1]
            w = bin[2]
            z = u + 1j * v
            self.u.append(u)
            self.v.append(v)
            self.w.append(w)
            self.z.append(z)

            self.current_speed.append(abs(z))
            direction = math.atan2(z.real, z.imag)*180/math.pi
            self.current_direction.append(direction)

    def __calc_bottom_stats(self, tolerance=30):
        bottom_bin = 1
        intensity = self.data['echo_intensity']['data']
        for bin_prev, bin_curr in zip(intensity, intensity[1:]):
            bin_flag_count = 0
            for beam_prev, beam_curr in zip(bin_prev, bin_curr):
                bin_diff = abs(beam_curr - beam_prev)
                if bin_diff > tolerance:
                    bin_flag_count += 1

            if bin_flag_count < 2:
                bottom_bin += 1
            else:
                break

        bottom_stats = {}
        bottom_stats['bottom_bin'] = bottom_bin
        bin_1_distance = (
            (self.data['fixed_leader']['bin_1_distance'] + self.transducer_depth) / 100  # NOQA
        )
        bottom_stats['range_to_bottom'] = (
            bottom_stats['bottom_bin'] *
            self.data['fixed_leader']['depth_cell_length'] / 100.0 +
            bin_1_distance
        )
        bottom_stats['side_lobe_start'] = (
            int(
                math.cos(
                    self.data['fixed_leader']['beam_angle'] *
                    math.pi/180.0
                ) *
                bottom_stats['range_to_bottom']
            )
        )
        bottom_stats['last_good_bin'] = bottom_stats['side_lobe_start'] - 1
        bottom_stats['last_good_counter'] = int(
            bottom_stats['last_good_bin'] - 1
        )
        self.bottom_stats = bottom_stats
        self.last_good_bin = bottom_stats['last_good_bin']
        self.last_good_counter = bottom_stats['last_good_counter']

    def battery_flag(self):
        """
        QARTOD Test #1 Strongly Recommended
        For TRDI ADCP's the battery flag test can not be performed.
        Confirmation that battery information is not sent out in
        data files was recieved by e-mail on 01/24/2014 from
        TRDI Field Service (Wilbur Rotoni)
        """
        return battery_flag_test(self.data)

    def checksum_flag(self):
        """
        QARTOD Test #2 Required
        The fact that this code is running means that the checksum was tested
        and passed the test. If the checksum test fails, we will
        not load anything into the database.
        """
        return checksum_test(self.data)

    def bit_flag(self):
        """
        Not an official QARTOD test.  Checks special TRDI bit flag.
        """
        bit_flags = bit_test(self.data['variable_leader']['bit_result'])

        return bit_flags

    def orientation_flags(self, max_pitch=20, max_roll=20):
        """
        QARTOD Test #3 Required: orientation (pitch and roll) tests
        """
        pitch = self.data['variable_leader']['pitch']/100.0
        roll = self.data['variable_leader']['roll']/100.0
        return orientation_test(pitch, roll,
                                max_pitch, max_roll)

    def sound_speed_flags(self, sound_speed_min=1400, sound_speed_max=1600):
        """
        QARTOD Test 4 Required: Sound speed test
        """
        return (
            sound_speed_test(self.data['variable_leader']['speed_of_sound'],
                             sound_speed_min, sound_speed_max)
        )

    def noise_floor_flag(self):
        """
        QARTOD Test #5 Strongly Recommended
        For TRDI ADCP's noise floor data is not available.
        """
        return noise_floor_test(self.data)

    def signal_strength_flag(self):
        """
        QARTOD Test #6 Strongly Recommended
        For TRDI ADCP's noise floor data is not available.
        """
        return signal_strength_test(self.data)

    def signal_to_noise_flag(self):
        """
        QARTOD Test #7 Strongly Recommended
        For TRDI ADCP's noise floor data is not available.
        """
        return signal_to_noise_test(self.data)

    def correlation_magnitude_flags(self,
                                    good_tolerance=115,
                                    questionable_tolerance=64):
        """
        QARTOD Test #8 Strongly Recommended
        correlation magnitude test
        """
        correlation = self.data['correlation']['data'][:self.last_good_bin]  # NOQA
        return (
            correlation_magnitude_test(correlation,
                                       good_tolerance,
                                       questionable_tolerance)
        )

    def percent_good_flags(self, percent_good=21, percent_bad=17):
        """
        QARTOD Test #9 Required
        percent good test
        default limits on this test are > 21%, good: < 17%, bad
        limits derived from TRDI Spreadsheed based on our instruments and setup
        """
        one_bad_percent = imap(itemgetter(2),
                               self.data['percent_good']['data'][:self.last_good_bin])  # NOQA
        all_good_percent = imap(itemgetter(3),
                                self.data['percent_good']['data'][:self.last_good_bin])  # NOQA
        return (
            percent_good_test(one_bad_percent, all_good_percent,
                              percent_good, percent_bad)
        )

    def current_speed_flags(self, max_speed=150):
        """
        QARTOD Test #10 Required
        current speed test
        150 cm/s is the West Florida Shelf limit.  Adjust as necessary
        """

        return current_speed_test(self.current_speed[:self.last_good_counter])

    def current_direction_flags(self):
        """
        QARTOD Test #11 Required
        current direction test
        Verifies that the current direction is less than 360 degrees
        Negative values are made positive by adding 360 to them
        """

        return current_direction_test(self.current_direction[:self.last_good_counter])  # NOQA

    def horizontal_velocity_flags(self,
                                  max_u_vel=150, max_v_vel=150):
        """
        QARTOD Test #12 Required
        horizontal velocity test
        If EITHER u or v greater than 150 cm/s, velocity is bad
        150 cm/s is a local WFS limit
        """

        return horizontal_velocity_test(self.u[:self.last_good_counter],
                                        self.v[:self.last_good_counter],
                                        max_u_vel, max_v_vel)

    def vertical_velocity_flags(self, max_w_velocity=15):
        """
        QARTOD Test #13 Strongly Recommended
        vertical velocity test
        if w greater than 15 cm/s (10% MAX speed from TRDI), w is bad
        """

        return vertical_velocity_test(self.w[:self.last_good_counter],
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

        error_velocities = imap(itemgetter(3),
                                self.data['velocity']['data'][:self.last_good_counter])  # NOQA
        return error_velocity_test(error_velocities,
                                   questionable_error_velocity,
                                   bad_error_velocity)

    def stuck_sensor_flag(self):
        """
        QARTOD Test #15 Strongly Recommended
        This test will not be performed. For as long as I have been involved
        in QARTOD, real time tests were meant to be performed on the current
        sample as if it were the only sample. This test requires examinimg the
        historical samples and is therefore NOT a real time data test.
        """

        return stuck_sensor_test(self.data)

    def echo_intensity_flags(self, tolerance=2):
        """
        QARTOD Test #16 Required
        echo intensity test
        """

        echo_intensities = self.data['echo_intensity']['data'][:self.last_good_counter]  # NOQA
        return echo_intensity_test(echo_intensities)

    def range_drop_off_flags(self, drop_off_limit=60):
        """
        QARTOD Test #17 Strongly Recommended
        range drop-off test
        Range Limit set to 60 as recommended in QARTOD spreadsheet (CO-OPS).
        The QARTOD recommended cut-off is 30. ???
        """

        echo_intensities = self.data['echo_intensity']['data'][:self.last_good_counter]  # NOQA
        return range_drop_off_test(echo_intensities, drop_off_limit)

    def current_speed_gradient_flags(self, tolerance=6):
        """
        QARTOD Test #18 Strongly Recommended
        current speed gradient test
        """

        return current_speed_gradient_test(self.current_speed,
                                           tolerance)

import sys
import argparse
from trdi_adcp_readers.readers import read_PD0_file


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_path", help="Path of PD0 file to parse")
    args = parser.parse_args()

    pd0_data = read_PD0_file(args.input_path, 0)

    qaqc = TRDIQAQC(pd0_data)

    print qaqc.bottom_stats
    print qaqc.echo_intensity_flags()

    return 1

if __name__ == '__main__':
    sys.exit(main())
