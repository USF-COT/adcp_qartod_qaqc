# Tests.py - A library of ADCP Quality Assurance and Quality Control tests
#
# Most comments refer to:
# "Manual for Real-Time Quality Control of In-Situ Current Observations"
# Integrated Ocean Observing Systems
# Quality Assurance of Real-Time Oceanographic Data
#
# Depends on numpy
#
# By: Jeff Donovan <jdonovan@usf.edu>
#     Michael Lindemuth <mlindemu@usf.edu>


import numpy as np

ADCP_FLAGS = {'good': 1, 'questionable': 2, 'bad': 3}


def battery_flag_test(ensemble):
    """
    QARTOD Test #1 Strongly Recommended
    battery flag test can not be performed ... if the GOES header f_code=G and
    the data can be decoded, then sample is good

    TODO: IMPLEMENT
    """
    return ADCP_FLAGS['good']


def checksum_test(ensemble):
    """
    QARTOD Test #2 Required
    Checksum test can not be performed because data is converted to/from
    Psuedo-ASCII during the GOES transmission... if the GOES header
    f_code=G and the data can be decoded, then sample is good

    TODO: IMPLEMENT
    """
    return ADCP_FLAGS['good']


def bit_test(ensemble_VL):
    """
    Not an official QARTOD test.  Checks special TRDI bit flag.
    """
    bit = ensemble_VL[9]
    if bit == '0':
        return ADCP_FLAGS['good']
    else:
        return ADCP_FLAGS['bad']


def orientation_test(pitches, rolls, max_pitch=20, max_roll=20):
    """
    QARTOD Test #3 Required: orientation (pitch and roll) tests
    """
    tilt_flags = []
    for pitch, roll in zip(pitches, rolls):
        if (abs(pitch) < max_pitch and abs(roll) < max_roll):
            tilt_flags.append(ADCP_FLAGS['good'])
        else:
            tilt_flags.append(ADCP_FLAGS['bad'])

    return tilt_flags


def sound_speed_test(ensemble_VL, sound_speed_min=1400, sound_speed_max=1600):
    """
    QARTOD Test 4 Required: Sound speed test
    """
    ssv = ensemble_VL[10]
    if (ssv > sound_speed_min and ssv < sound_speed_max):
        return ADCP_FLAGS['good']
    else:
        return ADCP_FLAGS['bad']

# TODO: Implement Tests 5, 6, and 7 for other ADCP types
# NOTE: QARTOD Tests 5, 6, and 7 cannot be performed on TRDI ADCP


def correlation_magnitude_test(ensemble_correlation,
                               good_tolerance=115, questionable_tolerance=64):
    """
    QARTOD Test #8 Strongly Recommended
    correlation magnitude test
    """
    bin_flags = []
    for bin in ensemble_correlation:
        counts = {'good': 0, 'questionable': 0, 'bad': 0}
        for beam_correlation in bin:
            if beam_correlation >= good_tolerance:
                counts['good'] += 1
            elif beam_correlation >= questionable_tolerance:
                counts['questionable'] += 1
            else:
                counts['bad'] += 1

        if counts['good'] == len(bin):
            bin_flags.append(ADCP_FLAGS['good'])
        elif counts['good'] + counts['questionable'] >= 3:
            bin_flags.append(ADCP_FLAGS['questionable'])
        else:
            bin_flags.append(ADCP_FLAGS['bad'])

    return bin_flags


def percent_good_test(pg3, pg4, last_good_bin,
                      percent_good=21, percent_bad=17):
    """
    QARTOD Test #9 Required
    percent good test
    default limits on this test are > 21%, good: < 17%, bad
    limits derived from TRDI Spreadsheet based on our instruments and setup
    """
    pg_flag = (np.ones(last_good_bin)) * 2
    pg_sum = (pg3[:last_good_bin]
              + pg4[:last_good_bin])
    pg_good = np.where(pg_sum > percent_good)
    pg_bad = np.where(pg_sum < percent_bad)
    pg_flag[pg_good] = ADCP_FLAGS['good']
    pg_flag[pg_bad] = ADCP_FLAGS['bad']

    return pg_flag


def current_speed_test(current_speed, last_good_bin, max_speed=150):
    """
    QARTOD Test #10 Required
    current speed test
    150 cm/s is the West Florida Shelf limit.  Adjust as necessary
    """

    current_speed_flags = np.ones(last_good_bin)
    good_speed = current_speed[:last_good_bin]
    bad_speed_index = np.where(good_speed > max_speed)
    current_speed_flags[bad_speed_index] = 3

    return current_speed_flags


def current_direction_test(current_direction, last_good_bin):
    """
    QARTOD Test #11 Required
    current direction test
    Verifies that the current direction is less than 360 degrees
    Negative values are made positive by adding 360 to them
    """

    current_direction_flags = np.ones(last_good_bin)

    current_good_direction = current_direction[:last_good_bin]
    current_direction_negative = np.where(current_good_direction < 0.0)
    current_good_direction[current_direction_negative] += 360

    bad_direction_index = np.where(current_good_direction > 360)
    current_direction_flags[bad_direction_index] = 3

    return current_direction_flags


def horizontal_velocity_test(u, v, last_good_bin,
                             max_u_velocity=150, max_v_velocity=150):
    """
    QARTOD Test #12 Required
    horizontal velocity test
    If EITHER u or v greater than 150 cm/s, velocity is bad
    150 cm/s is a local WFS limit
    """

    hvel_flags = np.ones(last_good_bin)

    u_index = np.where(abs(u[:last_good_bin]) > max_u_velocity)
    v_index = np.where(abs(v[:last_good_bin]) > max_v_velocity)
    hvel_flags[u_index] = 3
    hvel_flags[v_index] = 3

    return hvel_flags


def vertical_velocity_test(w, last_good_bin,
                           max_w_velocity=15):
    """
    QARTOD Test #13 Strongly Recommended
    vertical velocity test
    if w greater than 15 cm/s (10% MAX speed from TRDI), w is bad
    """

    vvel_flags = np.ones(last_good_bin)
    w_index = np.where(abs(w[:last_good_bin]) > 15.0)
    vvel_flags[w_index] = 3

    return vvel_flags


def error_velocity_test(error_velocities, last_good_bin,
                        questionable_error_velocity=2.6,
                        bad_error_velocity=5.2):
    """
    QARTOD Test #14 Required
    error velocity test
    err < 2.6 cm/s, good: err > 5.2 cm/s, bad
    limits derrived from TRDI Spreadsheed based on our instruments and setup
    """

    errvel_flag = (np.ones(last_good_bin))*2
    err_good = np.where(abs(error_velocities[:last_good_bin])
                        < questionable_error_velocity)
    err_bad = np.where(abs(error_velocities[:last_good_bin])
                       < bad_error_velocity)
    errvel_flag[err_good] = 1
    errvel_flag[err_bad] = 3

    return errvel_flag


def echo_intensity_test(echo_intensities, tolerance=2):
    """
    QARTOD Test #15 Required
    echo intensity test
    """

    echo_int_flags = [1]
    for bin_prev, bin_curr in zip(echo_intensities, echo_intensities[1:]):
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

    return echo_int_flags


def range_drop_off_test(echo_intensities, drop_off_limit=60):
    """
    QARTOD Test #16 Strongly Recommended
    range drop-off test
    Range Limit set to 60 as recommended in QARTOD spreadsheet (CO-OPS).
    The QARTOD recommended cut-off is 30. ???
    """

    range_flags = []
    for bin in echo_intensities:
        bin_flag_count = 0
        for beam in bin:
            if beam < drop_off_limit:
                bin_flag_count += 1

        if bin_flag_count >= 2:
            range_flags.append(3)
        else:
            range_flags.append(1)

    return range_flags


def current_speed_gradient_test(current_speed, last_good_bin,
                                tolerance=6):
    """
    QARTOD Test #17 Strongly Recommended
    current speed gradient test
    """

    speed_flags = [1]
    for speed_prev, speed_curr in zip(current_speed[:last_good_bin],
                                      current_speed[1:last_good_bin]):
        speed_diff = abs(speed_curr - speed_prev)
        if speed_diff <= tolerance:
            speed_flags.append(1)
        else:
            speed_flags.append(3)

    return speed_flags
