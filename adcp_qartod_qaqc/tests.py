# Tests.py - A library of ADCP Quality Assurance and Quality Control tests
#
# Most comments refer to:
# "Manual for Real-Time Quality Control of In-Situ Current Observations"
# Integrated Ocean Observing Systems
# Quality Assurance of Real-Time Oceanographic Data
#
# By: Jeff Donovan <jdonovan@usf.edu>
#     Michael Lindemuth <mlindemu@usf.edu>

from itertools import izip


ADCP_FLAGS = {
    'good': 1,
    'no_test': 2,
    'suspect': 3,
    'bad': 4,
    'missing_data': 9
}


def battery_flag_test(ensemble):
    """
    QARTOD Test #1 Strongly Recommended
    For TRDI ADCP's the battery flag test can not be performed. Confirmation
    that battery information is not sent out in data files was recieved by
    e-mail on 01/24/2014 from TRDI Field Service (Wilbur Rotoni)
    """

    return ADCP_FLAGS['no_test']


def checksum_test(ensemble):
    """
    QARTOD Test #2 Required
    The fact that this code is running means that the checksum was tested
    and passed the test. If the checksum test fails, we will not load anything
    into the database.
    """

    return ADCP_FLAGS['good']


def bit_test(bit_flag):
    """
    Not an official QARTOD test.  Checks special TRDI bit flag.
    """

    if bit_flag == '0':
        return ADCP_FLAGS['good']
    else:
        return ADCP_FLAGS['bad']


def orientation_test(pitch, roll, max_pitch=20, max_roll=20):
    """
    QARTOD Test #3 Required: orientation (pitch and roll) tests
    """

    if (abs(pitch) < max_pitch and abs(roll) < max_roll):
        tilt_flag = ADCP_FLAGS['good']
    else:
        tilt_flag = ADCP_FLAGS['bad']

    return tilt_flag


def sound_speed_test(sound_speed_velocity,
                     sound_speed_min=1400, sound_speed_max=1600):
    """
    QARTOD Test 4 Required: Sound speed test
    """

    if (sound_speed_velocity >= sound_speed_min and sound_speed_velocity <= sound_speed_max):  # NOQA
        return ADCP_FLAGS['good']
    else:
        return ADCP_FLAGS['bad']


def noise_floor_test(ensemble):
    """
    QARTOD Test #5 Strongly Recommended
    For TRDI ADCP's noise floor data is not available.
    """

    return ADCP_FLAGS['no_test']


def signal_strength_flag(ensemble):
    """
    QARTOD Test #6 Strongly Recommended
    For TRDI ADCP's signal strength data is not available.
    """

    return ADCP_FLAGS['no_test']


def signal_to_noise_test(ensemble):
    """
    QARTOD Test #7 Strongly Recommended
    For TRDI ADCP's signal to noise data is not available.
    """

    return ADCP_FLAGS['no_test']


def correlation_magnitude_test(ensemble_correlation,
                               good_tolerance=115, suspect_tolerance=64):
    """
    QARTOD Test #8 Strongly Recommended
    correlation magnitude test
    limits derived from TRDI Spreadsheet based on our instruments and setup
    """

    bin_flags = []
    for bin in ensemble_correlation:
        counts = {'good': 0, 'suspect': 0, 'bad': 0}
        for beam_correlation in bin:
            if beam_correlation >= good_tolerance:
                counts['good'] += 1
            elif beam_correlation >= suspect_tolerance:
                counts['suspect'] += 1
            else:
                counts['bad'] += 1

        if counts['good'] == len(bin):
            bin_flags.append(ADCP_FLAGS['good'])
        elif counts['good'] + counts['suspect'] >= 3:
            bin_flags.append(ADCP_FLAGS['suspect'])
        else:
            bin_flags.append(ADCP_FLAGS['bad'])

    return bin_flags


def percent_good_test(one_bad_percent_data, all_good_percent_data,
                      percent_good=21, percent_bad=17):
    """
    QARTOD Test #9 Required
    percent good test
    default limits on this test are > 21%, good: < 17%, bad
    limits derived from TRDI Spreadsheet based on our instruments and setup
    """

    pg_flags = []
    for one_bad_percent, all_good_percent in izip(one_bad_percent_data,
                                                  all_good_percent_data):
        pg_sum = one_bad_percent + all_good_percent
        if pg_sum >= percent_good:
            pg_flags.append(ADCP_FLAGS['good'])
        elif pg_sum <= percent_bad:
            pg_flags.append(ADCP_FLAGS['bad'])
        else:
            pg_flags.append(ADCP_FLAGS['suspect'])

    return pg_flags


def current_speed_test(current_speed, max_speed=150):
    """
    QARTOD Test #10 Required
    current speed test
    150 cm/s is the West Florida Shelf limit.  Adjust as necessary
    """

    current_flags = []
    for speed in current_speed:
        if speed <= max_speed:
            current_flags.append(ADCP_FLAGS['good'])
        else:
            current_flags.append(ADCP_FLAGS['bad'])

    return current_flags


def current_direction_test(current_direction):
    """
    QARTOD Test #11 Required
    current direction test
    Verifies that the current direction is less than 360 degrees
    Negative values are made positive by adding 360 to them
    """

    flags = []
    for direction in current_direction:
        if direction < 0.0:
            direction += 360

        if direction <= 360:
            flags.append(ADCP_FLAGS['good'])
        else:
            flags.append(ADCP_FLAGS['bad'])

    return flags


def horizontal_velocity_test(u, v,
                             max_u_velocity=150, max_v_velocity=150):
    """
    QARTOD Test #12 Required
    horizontal velocity test
    If EITHER u or v greater than 150 cm/s, velocity is bad
    150 cm/s is a local WFS limit
    """

    flags = []
    for u_vel, v_vel in izip(u, v):
        if abs(u_vel) > max_u_velocity or abs(v_vel) > max_v_velocity:
            flags.append(ADCP_FLAGS['bad'])
        else:
            flags.append(ADCP_FLAGS['good'])

    return flags


def vertical_velocity_test(w, max_w_velocity=15):
    """
    QARTOD Test #13 Strongly Recommended
    vertical velocity test
    if w greater than 15 cm/s (10% MAX speed from WFS), w is bad
    """

    flags = []
    for vel in w:
        if abs(vel) <= max_w_velocity:
            flags.append(ADCP_FLAGS['good'])
        else:
            flags.append(ADCP_FLAGS['bad'])

    return flags


def error_velocity_test(error_velocities,
                        suspect_error_velocity=2.6,
                        bad_error_velocity=5.2):
    """
    QARTOD Test #14 Required
    error velocity test
    err < 2.6 cm/s, good: err > 5.2 cm/s, bad
    limits derrived from TRDI Spreadsheet based on our instruments and setup
    """

    flags = []
    for error_vel in error_velocities:
        if error_vel < suspect_error_velocity:
            flags.append(ADCP_FLAGS['good'])
        elif error_vel < bad_error_velocity:
            flags.append(ADCP_FLAGS['suspect'])
        else:
            flags.append(ADCP_FLAGS['bad'])

    return flags


def stuck_sensor_test(ensemble):
    """
    QARTOD Test #15 Strongly Recommended
    This test will not be performed. For as long as I have been involved
    in QARTOD, real time tests were meant to be performed on the current
    sample as if it were the only sample. This test requires examinimg the
    historical samples and is therefore NOT a real time data test.
    """

    return ADCP_FLAGS['no_test']


def echo_intensity_test(echo_intensities, tolerance=2):
    """
    QARTOD Test #16 Required
    echo intensity test
    """

    flags = [1]
    for bin_prev, bin_curr in zip(echo_intensities, echo_intensities[1:]):
        bin_flag_count = 0
        for beam_prev, beam_curr in zip(bin_prev, bin_curr):
            beam_diff = int(beam_prev) - int(beam_curr)
            if beam_diff < tolerance:
                bin_flag_count += 1

        if bin_flag_count == 0:
            flags.append(ADCP_FLAGS['good'])
        elif bin_flag_count == 1:
            flags.append(ADCP_FLAGS['suspect'])
        else:
            flags.append(ADCP_FLAGS['bad'])

    return flags


def range_drop_off_test(echo_intensities, drop_off_limit=60):
    """
    QARTOD Test #17 Strongly Recommended
    range drop-off test
    Range Limit set to 60 as recommended in QARTOD spreadsheet (CO-OPS).
    The QARTOD recommended cut-off is 30. ???
    """

    flags = []
    for bin in echo_intensities:
        bin_flag_count = 0
        for beam in bin:
            if beam < drop_off_limit:
                bin_flag_count += 1

        if bin_flag_count >= 2:
            flags.append(ADCP_FLAGS['bad'])
        else:
            flags.append(ADCP_FLAGS['good'])

    return flags


def current_speed_gradient_test(current_speed,
                                tolerance=6):
    """
    QARTOD Test #18 Strongly Recommended
    current speed gradient test
    """

    flags = [1]
    for speed_prev, speed_curr in zip(current_speed,
                                      current_speed[1:]):
        speed_diff = abs(speed_curr - speed_prev)
        if speed_diff <= tolerance:
            flags.append(ADCP_FLAGS['good'])
        else:
            flags.append(ADCP_FLAGS['bad'])

    return flags
