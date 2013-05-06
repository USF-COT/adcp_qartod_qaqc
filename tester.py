
from pycurrents.adcp.rdiraw import Multiread
import numpy as np
import numpy.ma as ma
from datetime import datetime

class TRDI_ADCP_QARTOD_Tester(object):
    """
    Performs QARTOD Data Quality Assurance and Control Tests on TRDI ADCP Data according to:
    "IOOS: Manual for Real-Time Quality Control of In-Situ Current Observations"

    Expects data in the same format as output by University of Hawaii's Multiread function

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
        return TRDI_ADCP_QARTOD_Tester(m.read(), transducer_depth)

    def __init__(self, multiread_data, transducer_depth):
        self.data = multiread_data
        self.transducer_depth = transducer_depth

        self.__read_timestamp()
        self.__read_configuration()

        self.find_ensemble_bottom_bins()

    def __read_timestamp(self):
        yr=str(self.data.yearbase)
        mon=str(self.data.VL[0][2]).zfill(2)
        day=str(self.data.VL[0][3]).zfill(2)
        hr=str(self.data.VL[0][4]).zfill(2)
        min=str(self.data.VL[0][5]).zfill(2)
        sec=str(self.data.VL[0][6]).zfill(2)
        dt='/'.join([yr, mon, day])
        tm=':'.join([hr, min, sec])
        time_str=' '.join([dt, tm])
        self.timestamp = datetime.strptime(time_str, '%Y/%m/%d %H:%M:%S')

    def __read_configuration(self):
        """
        Reads some shared configuration variables into more easily accessible class variables
        """

        self.NCells = self.data.NCells # number of data bins
        self.raw_depth = self.data.dep # vector containing the uncorrected depth of each bin
        self.bin1depth = self.data.Bin1Dist + self.transducer_depth # actual depth of bin 1
        self.depth = self.raw_depth + self.transducer_depth # vector containing the corrected depth of each bin
        self.frequency = self.data.sysconfig['kHz'] # TRDI ADCP operating frequency
        self.beam_angle = self.data.sysconfig['angle'] # beam angle is necessary for calculating last good bin
        if (self.data.sysconfig['up'] == True): # orientation is either up of down and stored as a string
            self.orientation = 'up'
        else:
            self.orientation = 'down'
        self.bin_size = self.data.FL.CellSize # this is in cm instead of m
        self.blank = self.data.FL.Blank # this is in cm instead of m
        self.pings = self.data.FL.NPings # number of pings per ensemble (sample)
        self.transmit_lag = self.data.FL.TransLag # lag is used to calculate distance to first good bin
        self.transmit_pulse = self.data.FL.Pulse # pulse length is used to calculate distance to first good bin
        self.TPMin = str(self.data.FL.TPP_min).zfill(2) # ping interval - minutes
        self.TPSec = str(self.data.FL.TPP_sec).zfill(2) # ping interval - seconds
        self.TP100 = str(self.data.FL.TPP_hun).zfill(2) # ping interval - hundreds
        self.ping_interval = ':'.join([self.TPMin, self.TPSec, self.TP100]) # ping interval text string
        self.mag_declination = self.data.FL.EV / 10.0 # magnetic declination applied to data, site specific, user input

    def __read_velocities(self):
        """
        extract velocity from the masked array
        Univ of Hawaii code converts the mm/s velocity data from the instrument to
        m/s. We work in cm/s so I will have to convert the velocities from m/s to cm/s
        """
        u_vel_masked=self.data.vel1
        u_vel_masked.mask=ma.nomask
        u_vel=u_vel_masked.compressed()
        self.u = u_vel * 100.

        v_vel_masked=self.data.vel2
        v_vel_masked.mask=ma.nomask
        v_vel=v_vel_masked.compressed()
        self.v = v_vel * 100.

        w_vel_masked=self.data.vel3
        w_vel_masked.mask=ma.nomask
        w_vel=w_vel_masked.compressed()
        self.w = w_vel * 100.

        err_vel_masked=self.data.vel4
        err_vel_masked.mask=ma.nomask
        err_vel=err_vel_masked.compressed()
        self.ev = err_vel * 100.

        self.z = self.u + 1j + self.v
        self.current_speed = abs(self.z)
        self.current_direction = np.arctan2(self.z.real, self.z.imag)*180/np.pi

    def find_ensemble_bottom_bins(self, tolerance=30):
        self.ensemble_bottoms = []
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
            bottom_stats['range_to_bottom'] = (bottom_stats['bottom_bin'] * (self.NCells / 100)) + self.data.Bin1Depth
            bottom_stats['side_lobe_start'] = int(np.cos(self.data.sysconfig['angle'] * (np.pi/180.)) 
                * bottom_stats['range_to_bottom'])
            bottom_stats['last_good_bin'] = bottom_stats['side_lobe_start'] - 1
            bottom_stats['last_good_counter'] = bottom_stats['last_good_bin'] - 1
            self.ensemble_bottom_stats.append(bottom_stats)

    def battery_flag_test(self):
        """
        QARTOD Test #1 Strongly Recommended
        battery flag test can not be performed ... if the GOES header f_code=G and the data can be
        decoded, then sample is good
        """
        return 1

    def checksum_test(self):
        """
        QARTOD Test #2 Required
        Checksum test can not be performed because data is converted to/from Psuedo-ASCII during the
        GOES transmission... if the GOES header f_code=G and the data can be decoded, then sample is
        good
        """
        return 1

    def bit_test(self):
        """
        Not an official QARTOD test.  Checks special TRDI bit flag.
        """
        bit_flags = []
        for ensemble in self.data.VL:
            bit = ensemble[9]
            if bit == '0':
                bit_flags.append(1)
            else:
                bit_flags.append(3)

        return bit_flags

    def orientation_test(self, max_pitch=20, max_roll=20):
        """
        QARTOD Test #3 Required: orientation (pitch and roll) tests
        """
        tilt_flags = []
        for pitch, roll in zip(self.data.pitch, self.data.roll):
            if (abs(pitch) < max_pitch and abs(roll) < max_roll):
                tilt_flags.append(1)
            else:
                tilt_flags.append(3)

        return tilt_flags

    def sound_speed_test(self, sound_speed_min=1400, sound_speed_max=1600):
        """
        QARTOD Test 4 Required: Sound speed test
        """
        ssval_flags = []
        for ensemble_ssv in self.data.VL:
            ssv = ensemble_ssv[10]
            if (ssv > sound_speed_min and ssv < sound_speed_max):
                ssval_flags.append(1)
            else:
                ssval_flags.append(3)

        return ssval_flags

    # NOTE: QARTOD Tests 5, 6, and 7 cannot be performed on TRDI ADCP

    def correlation_magnitude_test(self, good_tolerance=115, questionable_tolerance=64):
        """
        QARTOD Test #8 Strongly Recommended
        correlation magnitude test
        """
        ensemble_flags = []
        for ensemble in self.data.corr:
            bin_flags = []
            for bin in ensemble:
                counts = {'good': 0, 'questionable': 0, 'bad': 0}
                for beam_correlation in bin:
                    if beam_correlation >= good_tolerance:
                        counts['good'] += 1
                    elif beam_correlation >= questionable_tolerance:
                        counts['questionable'] += 1
                    else:
                        counts['bad'] += 1
                if counts['good'] == len(bin):
                    bin_flags.append(1)
                elif counts['good'] + counts['questionable'] >= 3:
                    bin_flags.append(2)
                else:
                    bin_flags.append(3)
            ensemble_flags.append(bin_flags)

        return ensemble_flags

    def percent_good_test(self, percent_good=21, percent_bad=17):
        """
        QARTOD Test #9 Required
        percent good test
        default limits on this test are > 21%, good: < 17%, bad
        limits derived from TRDI Spreadsheed based on our instruments and setup
        """
        ensemble_flags = []
        for pg3, pg4, bottom_stats in zip(self.data.pg3, self.data.pg4, self.ensemble_bottom_stats):
            pg_flag = (np.ones(bottom_stats['last_good_counter'])) * 2
            pg_sum = pg3[:bottom_stats['last_good_counter']] + pg4[:bottom_stats['last_good_counter']]
            pg_good = np.where(pg_sum > percent_good)
            pg_bad = np.where(pg_sum < percent_bad)
            pg_flag[pg_good] = 1
            pg_flag[pg_bad] = 3
            ensemble_flags.append(pg_flag)

        return ensemble_flags

    def current_speed_test(self, max_speed=150):
        """
        QARTOD Test #10 Required
        current speed test
        150 cm/s is the West Florida Shelf limit.  Adjust as necessary
        """
        
        last_good_counter = self.ensemble_bottom_stats[0]['last_good_counter']
        current_speed_flags = np.ones(last_good_counter)
        good_speed = self.current_speed[:last_good_counter]
        bad_speed_index = np.where(good_speed > max_speed)
        current_speed_flags[bad_speed_index] = 3

        return current_speed_flags

    def current_direction_test(self):
        """
        QARTOD Test #11 Required
        current direction test
        Verifies that the current direction is less than 360 degrees
        Negative values are made positive by adding 360 to them
        """

        last_good_counter = self.ensemble_bottom_stats[0]['last_good_counter']
        
        good_direction = self.current_direction[:last_good_counter]
        good_direction_negative = np.where(good_direction < 0.0)
        good_directon[current_direction_negative] += 360

        bad_direction_index = np.where(good_direction > 360)
        current_direction_flags[bad_direction_index] = 3
        
        return current_direction_flags

    def horizontal_velocity_test(self, max_u_velocity=150, mav_v_velocity=150):
        """
        QARTOD Test #12 Required
        horizontal velocity test
        If EITHER u or v greater than 150 cm/s, velocity is bad
        150 cm/s is a local WFS limit
        """
        
        last_good_counter = self.ensemble_bottom_stats[0]['last_good_counter']
        hvel_flags = np.ones(last_good_counter)

        u_index = np.where(abs(self.u[:last_good_counter]) > 150)
        v_index = np.where(abs(self.v[:last_good_counter]) > 150)
        hvel_flags[u_index] = 3
        hvel_flags[v_index] = 3

        return hvel_flags

    def vertical_velocity_test(self, max_w_velocity=15):
        """
        QARTOD Test #13 Strongly Recommended
        vertical velocity test
        if w greater than 15 cm/s (10% MAX speed from TRDI), w is bad
        """
        
        last_good_counter = self.ensemble_bottom_stats[0]['last_good_counter']
        vvel_flags = np.ones(last_good_counter)
        w_index = np.where(abs(self.w[:last_good_counter]) > 15.0)
        vvel_flags[w_index] = 3

        return vvel_flags

    def error_velocity_test(self, questionable_error_velocity=2.6, bad_error_velocity=5.2):
        """
        QARTOD Test #14 Required
        error velocity test
        err < 2.6 cm/s, good: err > 5.2 cm/s, bad
        limits derrived from TRDI Spreadsheed based on our instruments and setup
        """

        last_good_counter = self.ensemble_bottom_stats[0]['last_good_counter']
        errvel_flag = (np.ones(last_good_counter))*2
        err_good = np.where(abs(self.ev[:last_good_counter]) < questionable_error_velocity)
        err_bad = np.where(abs(self.ev[:last_good_coutner]) < bad_error_velocity)
        errvel_flag[err_good] = 1
        errvel_flag[err_bad] = 3

        return errvel_flag

    def echo_intensity_test(self, tolerance=2):
        """
        QARTOD Test #15 Required
        echo intensity test
        """

        ensemble_flags = []
        for ensemble in self.data.amp:
            echo_int_flags = [1]
            for bin_prev, bin_curr in zip(ensemble, ensemble[1:]):
                bin_flag_count = 0
                for beam_prev, beam_curr in zip(bin_prev, bin_curr):
                    beam_diff = int(beam_prev) - int(beam_curr)
                    if beam_diff < tolerance:
                        bin_flag_count += 1
                
                if bin_flag_count == 0:
                    echo_int_flags.append(1)
                elif bin_flag_count == 1:
                    echo_int_flags.append(2)
                else:
                    echo_int_flags.append(3)

            ensemble_flags.append(echo_int_flags)

        return ensemble_flags

    def range_drop_off_test(self, drop_off_limit=60):
        """
        QARTOD Test #16 Strongly Recommended
        range drop-off test
        Range Limit set to 60 as recommended in QARTOD spreadsheet (CO-OPS).
        The QARTOD recommended cut-off is 30. ???
        """

        ensemble_flags = []
        for ensemble in self.data.amp:
            range_flags = []
            for bin in ensemble:
                bin_flag_count = 0
                for beam in bin:
                    if beam < drop_off_limit:
                        bin_flag_count += 1

                if bin_flag_count >= 2:
                    range_flags.append(3)
                else:
                    range_flags.append(1)

            ensemble_flags.append(range_flags)

        return ensemble_flags

    def current_speed_gradient_test(self, tolerance=6):
        """
        QARTOD Test #17 Strongly Recommended
        current speed gradient test
        """
        
        last_good_counter = self.ensemble_bottom_stats[0]['last_good_counter']
        good_speed = self.current_speed[:last_good_counter]
        speed_flags = [1]
        for speed_prev, speed_curr in zip(good_speed, good_speed[1:]):
            speed_diff = abs(speed_curr - speed_prev)
            if speed_diff <= tolerance:
                speed_flags.append(1)
            else:
                speed_flags.append(3)

        return speed_flags

import sys
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_path", help="Path of PD0 file to parse")
    parser.add_argument("read_type", help="Multiread read type (e.g., wh).  See Mutiread documentation for details")
    parser.add_argument("transducer_height", type=float, help="Depth of ADCP transducer")
    args = parser.parse_args()

    trdi_qaqc = TRDI_ADCP_QARTOD_Tester.from_file(args.input_path, args.read_type, args.transducer_height)
    print 'Ensemble Bottom Bins: %s' % (trdi_qaqc.ensemble_bottoms,)
    
    return 1

if __name__ == '__main__':
    sys.exit(main())
