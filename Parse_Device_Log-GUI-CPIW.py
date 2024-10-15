#####################################################################
# Date Modified: 2023-07-12
# Version: 3.0.1.5 (CHANGE BELOW IN CODE AS WELL!!!!)
#
# This script will take and parse all ".dat", i.e. log files, in
# the folder where it is run and look for pegged acceleration data,
# Flat areas where Vector Magnitude is not close to 1g, as well as
# parse the defined record types.
#
# There is an option to write activity data to csv file(s).  If the
# log file has more than 28800 activity records (8 hours) it will
# break the csv files up and add the first timestamp in the file
# to the file name.
#
# There is also an option to only parse a certain date/time range
#
# Added option to calibrate by either the .cal from nebulizer OR
# .json from the agdc file.
########################################################################

PySimpleGUI_License = 'eeyXJAMPa9WvNWlabUnvNVljVNHSl4wAZXS6IT6YIBkSRfpEc43YRyyDadWCJG19djG4l3vybci4IDsjIVkIxHpzYd2CV1uKcy2mVJJtRRCiIB6tMwTRchytOgD9kE51OoDkAIxwOES2wIitTSGTlsj0ZwWn5yzbZ4UCRGlAc0G6xBv9epWO1mlJbenuRCWQZ2XZJVzCadW19XubI1jIo6iSNJSm4xwyI6iswSiyT0mnFbt8ZAU2ZCp8c5nwNd0bIXjXoDiMRzGf981qZ9yzIwsYIakF5Xh4bqWIV4MBYPX2Nn0GIqjcojiIQE3CJFvKcc3rM8iLLMCSJKD6b12f1twwYxWK5I50IGjrociHQwWyNh0FacUNddypYpXgBcoEIsiUwmi0QO3iVBzwdSGu9zt8ZBX0JbJXRzCvIa6BIejyQFxpNETxYg3zIeiYwuiuRbGgFd05ZRUclszrcU3gV1lSZrCVIv6HIuj7IuwRM9jPQ0t8MoTVAIt2MDT8UsiFLqCiJFEjYxXlRolBRNXYh0wjaeXuJqlQceySIR6nI6j7IpwvMdjVUAtPM4TOAJtBMLDDguiOLQCjJlFTb9WOF9pMbpETFFkHZGHoJclfcC37MKiJO0i1JGkxbl3JVynKLCm1NBycbH3gNdzVQMHvRCo2Z4W8FUjpdNGplQnKc0mFFrwHaGCD5djvbI2q0XiaLoCXJ7JiUhEjFPknZGHKJdlPcu3bM6iZOjiNIhzgOMCs45xSODD0Y6urNYjIYNu1MQTGMvi5fyQ7=L=o17e77197f7c0ced6424e953715a3f5ba1afa32b1f6a75652f3e48ade63e58598187eec86349cd58fa9c6d8d445e202b8bb2396401b65c0594d555c4e00024099acf5fe818588186849f854364e737e7781fedf082c322e1d84d2bb81ad76bba575c0968e7bd027603e0f90c3e5167d6e9bdf948faec9370977cf4497e3573d1b7e1f06ae609dd5c17db11c99d16bbaadcac71314117e4c8f824e6a2e93ff4f99c7549b0629b4f4ce59c70114a0045c13f45e0884281bc2d92a4d07ac24d169b61468c4c0b9db63701c4b79b12f7f07ea4c644f46e1e510a5e3135b563bb987fe46d36eed9ede9b09437b9a1e5d3d7c3715e9695013292617c4ac4f5fed5c2d9f803f6a47a2c2eb9bbf631ca188d174a8e1e30d266fcee0e91a21e396a51daa6d868c875a4ccc5d36fb1c531be0252f283e14f315fa1ee98b99f3bfbe9b73f72d749da33738674ecd2d1f38180b9da32f9db6c422a4b71677684435411d6e0bbccab7dd841dfdca0e942c3da1412152dd022a2e581b4187d28c241b8dcbdecfbba958c577134bf05385b67abf802b76661cc7ccf2c7295786a56dcf312f339aee3ff35ad5e1473366fcd956ee831b176b28d5a1954c1649a9206c56668bca3f112e36fcbc5eebd1fdc53ba019c0e9155801ce16b8bbde0417d076db456fff116e81c9229f9d16370bf92a3751ed37abca1e1fc64220fd61a45d150adfd9749c06'

from pathlib import Path

import PySimpleGUI as sg
import os.path
import fnmatch
from logparser3_9 import parse
import os
from datetime import datetime, timezone, timedelta
from bitstring import BitStream
import math
import numpy as np
import csv
from threading import Thread
import json
import pandas as pd          # library to get function to read csv into a data frame
import plotly.express as px  # this library is to do plotting
import plotly.io as pio      # this is to direct where plot output goes
import zipfile
import struct

pio.renderers.default = 'browser'  # this is to plot into default web broswer


VERSION = '3.0.1.5'   # CHANGE ABOVE IN FILE DESCRIPTION!!!

# ###########
# # OPTIONS #
# ###########
# LOG_ACTIVITY_DATA = True
# NO_FILTER = 0
# USE_CAL_VALUES = False


##########################################################################

#####################
# CONSTANTS/GLOBALS #
#####################
message = ''
progress = 0
totallogsize = 0

File_Split_Level = 172800  # Seconds
consecutiveSamples = 8
LOG_IMU = 0
createHtmlPlot = 0

# Filter Coefficient for Flat area detection #
alpha = 0.8

# Default Sample Rate #
# FS = 32

# Date/Time Format #
FMT = '%Y/%m/%d %H:%M:%S'
CalibrationOrder = ["negativeZeroGOffsetX_32",
                    "negativeZeroGOffsetY_32",
                    "negativeZeroGOffsetZ_32",
                    "positiveZeroGOffsetX_32",
                    "positiveZeroGOffsetY_32",
                    "positiveZeroGOffsetZ_32",
                    "zeroGOffsetX_32",
                    "zeroGOffsetY_32",
                    "zeroGOffsetZ_32",
                    "offsetX_32",
                    "offsetY_32",
                    "offsetZ_32",
                    "sensitivityXX_32",
                    "sensitivityYY_32",
                    "sensitivityZZ_32",
                    "sensitivityXY_32",
                    "sensitivityXZ_32",
                    "sensitivityYZ_32",
                    "negativeZeroGOffsetX_64",
                    "negativeZeroGOffsetY_64",
                    "negativeZeroGOffsetZ_64",
                    "positiveZeroGOffsetX_64",
                    "positiveZeroGOffsetY_64",
                    "positiveZeroGOffsetZ_64",
                    "zeroGOffsetX_64",
                    "zeroGOffsetY_64",
                    "zeroGOffsetZ_64",
                    "offsetX_64",
                    "offsetY_64",
                    "offsetZ_64",
                    "sensitivityXX_64",
                    "sensitivityYY_64",
                    "sensitivityZZ_64",
                    "sensitivityXY_64",
                    "sensitivityXZ_64",
                    "sensitivityYZ_64",
                    "negativeZeroGOffsetX_128",
                    "negativeZeroGOffsetY_128",
                    "negativeZeroGOffsetZ_128",
                    "positiveZeroGOffsetX_128",
                    "positiveZeroGOffsetY_128",
                    "positiveZeroGOffsetZ_128",
                    "zeroGOffsetX_128",
                    "zeroGOffsetY_128",
                    "zeroGOffsetZ_128",
                    "offsetX_128",
                    "offsetY_128",
                    "offsetZ_128",
                    "sensitivityXX_128",
                    "sensitivityYY_128",
                    "sensitivityZZ_128",
                    "sensitivityXY_128",
                    "sensitivityXZ_128",
                    "sensitivityYZ_128",
                    "negativeZeroGOffsetX_256",
                    "negativeZeroGOffsetY_256",
                    "negativeZeroGOffsetZ_256",
                    "positiveZeroGOffsetX_256",
                    "positiveZeroGOffsetY_256",
                    "positiveZeroGOffsetZ_256",
                    "zeroGOffsetX_256",
                    "zeroGOffsetY_256",
                    "zeroGOffsetZ_256",
                    "offsetX_256",
                    "offsetY_256",
                    "offsetZ_256",
                    "sensitivityXX_256",
                    "sensitivityYY_256",
                    "sensitivityZZ_256",
                    "sensitivityXY_256",
                    "sensitivityXZ_256",
                    "sensitivityYZ_256"
                    ]

TemperatureCalibrationOrder = [ "mcuTempHigh",
                                "mcuTempLow",
                                "adxlTempHigh",
                                "adxlTempLow",
                                "tempHigh",
                                "tempLow",
                                "mcuTempCal1",
                                "mcuTempCal2",
                                "calibrationMethod",
                                "calibrationTime",
                                "isCalibrated"
                                ]

#############
# FUNCTIONS #
#############


def check_for_peg(xp, yp, zp, activityCountp, calibratedp):
    if activityCountp <= 1:
        check_for_peg.xcount = check_for_peg.ycount = check_for_peg.zcount = 0
        check_for_peg.status = 0
    if calibratedp:
        maxvalue = 7.0
    else:
        maxvalue = 2000
    # Check For Pegged Data #
    if abs(xp) > maxvalue:
        check_for_peg.xcount = check_for_peg.xcount + 1
    elif abs(yp) > maxvalue:
        check_for_peg.ycount = check_for_peg.ycount + 1
    elif abs(zp) > maxvalue:
        check_for_peg.zcount = check_for_peg.zcount + 1
    else:
        check_for_peg.xcount = check_for_peg.ycount = check_for_peg.zcount = 0
        check_for_peg.status = 0

    if (
        check_for_peg.xcount > consecutiveSamples
        or check_for_peg.ycount > consecutiveSamples
        or check_for_peg.zcount > consecutiveSamples
    ):
        check_for_peg.status = 1

    return check_for_peg.status


check_for_peg.xcount = check_for_peg.ycount = check_for_peg.zcount = 0
check_for_peg.status = 0


def check_for_flat(xf, yf, zf, activityCountf, calibratedf, flatFS):
    if activityCountf <= 1:
        check_for_flat.y_xn = check_for_flat.y_yn = check_for_flat.y_zn = 0
        check_for_flat.previous_x = (
            check_for_flat.previous_y
        ) = check_for_flat.previous_z = 0
        check_for_flat.flatCount_x = (
            check_for_flat.flatCount_y
        ) = check_for_flat.flatCount_z = 0
        check_for_flat.flatStatus = check_for_flat.vmCount = 0

    if calibratedf:
        flatLimit = 0.05
        VectorMag = math.sqrt(
            check_for_flat.y_xn**2
            + check_for_flat.y_yn**2
            + check_for_flat.y_zn**2
        )
    else:
        flatLimit = 5
        VectorMag = (
            math.sqrt(
                check_for_flat.y_xn**2
                + check_for_flat.y_yn**2
                + check_for_flat.y_zn**2
            )
            / 250
        )

    check_for_flat.y_xn = alpha * xf + (1 - alpha) * check_for_flat.previous_x
    check_for_flat.y_yn = alpha * yf + (1 - alpha) * check_for_flat.previous_y
    check_for_flat.y_zn = alpha * zf + (1 - alpha) * check_for_flat.previous_z
    check_for_flat.y_vm = alpha * VectorMag + (1 - alpha) * check_for_flat.previous_vm

    if activityCountf > 1:
        # CHECK X #
        if abs(check_for_flat.y_xn - xf) < flatLimit:
            check_for_flat.flatCount_x += 1
        else:
            check_for_flat.flatCount_x = 0

        # CHECK Y #
        if abs(check_for_flat.y_yn - yf) < flatLimit:
            check_for_flat.flatCount_y += 1
        else:
            check_for_flat.flatCount_y = 0

        # CHECK Z #
        if abs(check_for_flat.y_zn - zf) < flatLimit:
            check_for_flat.flatCount_z += 1
        else:
            check_for_flat.flatCount_z = 0

        # CHECK VECTOR MAG #
        if check_for_flat.y_vm > 1.5:
            check_for_flat.vmCount += 1
        else:
            check_for_flat.vmCount = 0

        if (
            check_for_flat.flatCount_x == 0
            and check_for_flat.flatCount_y == 0
            and check_for_flat.flatCount_z == 0
        ):
            check_for_flat.flatStatus = 0

        if ((check_for_flat.flatCount_x > flatFS and check_for_flat.flatCount_y > flatFS)
            or (check_for_flat.flatCount_z > flatFS and check_for_flat.flatCount_y > flatFS)
            or (check_for_flat.flatCount_x > flatFS and check_for_flat.flatCount_z > flatFS)) \
                and check_for_flat.vmCount > flatFS:
            check_for_flat.flatStatus = 1

    check_for_flat.previous_x = xf
    check_for_flat.previous_y = yf
    check_for_flat.previous_z = zf
    check_for_flat.previous_vm = VectorMag

    return check_for_flat.flatStatus


check_for_flat.y_xn = (
    check_for_flat.y_yn
) = check_for_flat.y_zn = check_for_flat.y_vm = 0
check_for_flat.previous_x = (
    check_for_flat.previous_y
) = check_for_flat.previous_z = check_for_flat.previous_vm = 0
check_for_flat.flatCount_x = check_for_flat.flatCount_y = check_for_flat.flatCount_z = 0
check_for_flat.flatStatus = 0
check_for_flat.vmCount = 0


def import_calibration_values(cal_file):

    with open(cal_file, 'r') as cfile:
        v = 0
        if cal_file.lower().endswith("cal"):
            csvReader = csv.reader(cfile)
            for row in csvReader:
                import_calibration_values.calval_temp[v] = row[1]
                v += 1
        else:
            jsonReader = json.load(cfile)
            while v < len(CalibrationOrder):
                import_calibration_values.calval_temp[v] = jsonReader[
                    CalibrationOrder[v]
                ]
                v += 1
    return import_calibration_values.calval_temp


import_calibration_values.calval_temp = np.empty(72, int)

def import_temperature_calibration_values(temp_cal_file):

    with open(temp_cal_file, 'r') as cfile:
        v = 0
        jsonReader = json.load(cfile)
        while v < len(TemperatureCalibrationOrder):
            import_temperature_calibration_values.calval_temperature[v] = jsonReader[TemperatureCalibrationOrder[v]]
            v += 1
    return import_temperature_calibration_values.calval_temperature


import_temperature_calibration_values.calval_temperature = np.empty(12, int)

def ticks_to_unix(ticks):
    # Calculate the difference between Unix epoch and .NET epoch
    epoch_diff = datetime(1970, 1, 1) - datetime(1, 1, 1)
    epoch_diff_ticks = epoch_diff.total_seconds() * 10 ** 7

    # Calculate Unix timestamp
    ticks_to_unix_timestamp = (ticks - epoch_diff_ticks) / 10 ** 7
    return ticks_to_unix_timestamp

def main_process(Log_Activity_Data, NoFilter, UseCalValues, zipfiles, basePath, begin_timestamp, end_timestamp):

    ####### Create Folder Structure ########################################################
    adxl_dir = '/Activity Files - Primary Accel/'
    timeGap_dir = '/Time Gap Files/'
    battery_dir = '/Battery Logs/'
    temperature_dir = '/Temperature_Files/'
    calibration_dir = '/Calibration Logs/'
    output_dir = '/output_files/'
    time_stamp_only_dir = 'ts_only/'
    epoch_dir = '/Epoch Files/'

    outPutPath = basePath + output_dir

    Path(outPutPath + adxl_dir + time_stamp_only_dir).mkdir(parents=True, exist_ok=True)
    Path(outPutPath + timeGap_dir).mkdir(parents=False, exist_ok=True)
    Path(outPutPath + battery_dir).mkdir(parents=False, exist_ok=True)
    Path(outPutPath + calibration_dir).mkdir(parents=False, exist_ok=True)
    Path(outPutPath + temperature_dir).mkdir(parents=False, exist_ok=True)
    Path(outPutPath + epoch_dir).mkdir(parents=False, exist_ok=True)


    ###########################################################################################

    filenumber = 1
    temperatureCalibration = False

    global message, progress, totallogsize

    if begin_timestamp != "":
        begin_timestamp = datetime.fromtimestamp(int(begin_timestamp), tz=timezone.utc)

    if end_timestamp != "":
        end_timestamp = datetime.fromtimestamp(int(end_timestamp), tz=timezone.utc)

    d = datetime.now()
    with open(outPutPath + '/Parse_Summary_%s.txt' % d.strftime('%Y_%m_%d-%H_%M_%S'), "w") as sumFile:

        sumFile.write('APP_VERSION: %s\n\n' % VERSION)
        print('\nAPP_VERSION: %s\n' % VERSION)
        # Create and open Summary File #
        fout_summary = open(outPutPath + 'summary_file.csv', "w")
        fout_summary.write('filename,firmware version,First Activity Timestamp,Last Activity Timestamp,'
                           'Activity records,Unexpected Resets,Expected Resets,Pegs,flat areas,'
                           'Timestamp Gaps\n')

        #zipfiles = 1
        if zipfiles:
            files = Path(basePath).glob('*.[ag][gt][d3][cx]')
            filetotal = len(fnmatch.filter(os.listdir(basePath), '*.agdc'))
            filetotal += len(fnmatch.filter(os.listdir(basePath), '*.gt3x'))
        else:
            files = Path(basePath).glob('*.[db][ai][tn]')
            filetotal = len(fnmatch.filter(os.listdir(basePath), '*.bin'))
            filetotal += len(fnmatch.filter(os.listdir(basePath), '*.dat'))

        # Get all log files #
        #for file in glob.glob(basePath + "/*.[dat][bin]"):
        for file in files:
            skip_file = 0
            print("file %s of %s" % (filenumber, filetotal))
            filenumber += 1
            if file.name.startswith('CPW'):
                min_dateTime_unix = 1514764800
            else:
                min_dateTime_unix = 1262304000

            if zipfiles:
                # Iterate through each file #
                if zipfile.is_zipfile(file):
                    with zipfile.ZipFile(file, 'r') as zf:
                        zf.extract('log.bin', basePath)
                        if str(file).endswith('.agdc'):
                            zf.extract('calibration.json', basePath)
                            try:
                                zf.extract('temperature_calibration.json', basePath)
                                tempcalvals = import_temperature_calibration_values(basePath + '/temperature_calibration.json')
                                if tempcalvals[0] != 0:
                                    temperatureCalibration = True
                                    gainMcu = (tempcalvals[4] - tempcalvals[5]) / (tempcalvals[0] - tempcalvals[1])
                                    gainAdxl = (tempcalvals[4] - tempcalvals[5]) / (tempcalvals[2] - tempcalvals[3])
                                    offsetMcu = tempcalvals[4] - (gainMcu * tempcalvals[0])
                                    offsetAdxl = tempcalvals[4] - (gainAdxl * tempcalvals[2])
                                else:
                                    temperatureCalibration = False
                                    print('\nWARNING!!!!!!!   Temperature Calibration File has default values!!!!!!!!!\n')
                            except:
                                temperatureCalibration = False


                            zf.extract('info.json', basePath)

                            with open(basePath + '/info.json', 'r') as info_file:
                                jsonReader = json.load(info_file)
                                firmware_version = jsonReader["firmware"]
                                #target_start_time_unix = jsonReader["startDate"]
                                downloadDate_unix = jsonReader["lastSampleTime"] + 86400
                        else:
                            zf.extract('info.txt',basePath)

                            with open(basePath + '/info.txt', 'r') as info_file:
                                for line in info_file:
                                    if "Firmware" in line:
                                        firmware_version = line[-6:].rstrip('\n')
                                    if "Last Sample Time" in line:
                                        number = int(line[-19:].rstrip('\n'))
                                        downloadDate_unix = ticks_to_unix(number)

                    filetoopen = basePath + '/log.bin'
                else:
                    skip_file = 1
            else:
                filetoopen = file
                firmware_version = "dat or bin file only"

            if not skip_file:
                filename = Path(file).stem
                with open(filetoopen, 'rb') as fin:
                    if Log_Activity_Data:
                        fout = open(outPutPath + adxl_dir + filename + '.csv', "w")
                        fout.write('ts,t,x,y,z,vm\n')

                    adxl = open(outPutPath + adxl_dir + time_stamp_only_dir + filename + '_ts_only.csv', "w")
                    adxl.write('ts,record_length,Unix Timestamp\n')

                    # Create and open Timestamp Issue CSV File #
                    fout1 = open(outPutPath + timeGap_dir + filename
                                 + '_datetime_Gap.csv', "w")

                    # Create and open Battery Log CSV File #
                    fout2 = open(outPutPath + battery_dir + filename
                                 + 'battery_log.csv', "w")
                    fout2.write("Time Stamp,Batter_Voltage\n")

                    # Create and open Temperature Log CSV File #
                    fout3 = open(outPutPath + temperature_dir + filename
                                 + 'temperature_log.csv', "w")
                    fout3.write("Time Stamp,ADXL_Temp,STM32_Temp, STM32_CAL1, STM32CAL2\n")

                    # Create and open Epoch File Log CSV File #
                    fout4 = open(outPutPath + epoch_dir + filename
                                 + 'epoch.csv', "w")
                    fout4.write("Time Stamp,X, Y, Z\n")

                    # Create and open calibration Log CSV File #
                    fout_cal = open(outPutPath + calibration_dir + filename
                                 + 'calibration_log.csv', "w")
                    fout_cal.write('ts,x,y,z\n')

                    # Print File to console and Summary.txt file #
                    print('############  ' + file.name + '  ############')
                    if zipfiles:
                        print('firmware version: %s' % firmware_version)
                    sumFile.write('Filename: %s\n' % str(file))
                    if zipfiles:
                        sumFile.write('Firmware Version: %s\n' % firmware_version)
                    parse_start_date = datetime.now()
                    sumFile.flush()
                    os.fsync(sumFile)

                    # Initialize Variables at beginning of each file #
                    i = 0
                    activity_count = 0
                    expected_reset_count = 0
                    unexpected_reset_count = 0
                    peggedState = False
                    pegs = 0
                    flatRegionState = False
                    flats = 0
                    firstFile = True
                    lasttimestamp = datetime.fromtimestamp(94694400, tz=timezone.utc)  #January 1, 1973 12:00:00 AM
                    timestampGap = 0
                    isCalVals = False
                    isCalibrated = False
                    record_type = 255
                    isIMUdata = 0
                    timedelta = 0.0
                    timedelta_accel = 0.0
                    timedelta_temp = 0.0
                    incrementsize = 1
                    totallogsize = os.path.getsize(filetoopen)
                    cal_orienatation = False
                    numberRows = 0
                    invalidRecordCount = 0
                    firstTimestampFound = False
                    getInvalidRecordBytes = False

                    # Create Raw Activity CSV file if Logging is enabled #
                    if Log_Activity_Data:
                        fout1.write("Current Time Stamp,Previous Time Stamp, Delta\n")

                    ##########################
                    # Start Parsing Log file #
                    ##########################
                    if not(zipfiles):
                        downloadDate_unix = int((datetime.now(tz=timezone.utc)- datetime.fromtimestamp(0, timezone.utc)).total_seconds())

                    for record in parse(fin, min_dateTime_unix, downloadDate_unix):
                        timestamp = record.timestamp.strftime(FMT)
                        unixTime = (record.timestamp - datetime.fromtimestamp(0, timezone.utc)).total_seconds()

                        #progress += total / (update * 10)
                        if not (NoFilter) and record.timestamp > end_timestamp:
                            break

                        if NoFilter or (begin_timestamp < record.timestamp < end_timestamp):
                            ############################
                            # Invalid Record with      #
                            # Valid Timestamp          #
                            ############################
                            error_list = [253, 254]
                            if record.type in error_list or getInvalidRecordBytes is True:

                                if getInvalidRecordBytes is True:
                                    sumFile.write('%-35s %-10s\n' % ('Invalid Record Size: ', str(record.bad_size)))
                                    getInvalidRecordBytes = False
                                else:
                                    invalidRecordCount = invalidRecordCount + 1
                                    address = hex(record.payload)
                                    if record.type == 253:
                                        sumFile.write(
                                            '%-35s %-10s %10s\n' % ('Checksum Collision: ', timestamp,
                                                                    str(unixTime)))
                                        sumFile.write(
                                            '%-35s %-10s\n' % ('Checksum Collision @ Position: ', str(address)))
                                    else:
                                        sumFile.write(
                                            '%-35s %-10s %10s\n' % ('Invalid Record @ Timestamp: ', timestamp,
                                                                    str(unixTime)))
                                        sumFile.write(
                                            '%-35s %-10s\n' % ('Invalid Record @ Position: ', str(address)))
                                    getInvalidRecordBytes = True

                            eventList = [0, 26, 27, 37, 36, 35, 33, 32, 31, 100]

                            if record.type in eventList and firstTimestampFound is False:
                                firstTimestamp = record.timestamp
                                firstTimestampUnix = unixTime
                                print('%-35s %-10s %10s' % (
                                    'First Record Timestamp: ', timestamp, str(unixTime)))
                                sumFile.write(
                                    '%-35s %-10s %10s\n' % ('First Record Timestamp: ', timestamp,
                                                            str(unixTime)))
                                firstTimestampFound = True

                            ############################
                            # Parse if Activity Record #
                            ############################
                            if (record.type == 27 or record.type == 26 or record.type == 0):
                                record_type = record.type
                                activityRecordLen = len(record.payload) / 2

                                # Locate and log first Activity Timestamp #
                                if i < 1:
                                    if record.type == 27:
                                        base_file = os.path.splitext(file)[0]
                                        if os.path.isfile(base_file + '.cal') or os.path.isfile(basePath + '/calibration.json') and UseCalValues:
                                            FS = int(len(record.payload) / 2 / 4.5)
                                            isCalVals = True
                                            isCalibrated = True
                                            if os.path.isfile(base_file + '.cal'):
                                                calvals = import_calibration_values(base_file + '.cal')
                                            else:
                                                calvals = import_calibration_values(basePath + '/calibration.json')
                                            txt032 = calvals[0:18]
                                            txt064 = calvals[18:36]
                                            txt128 = calvals[36:54]
                                            txt256 = calvals[54:72]

                                            if FS == 32:
                                                txt = txt032
                                                O = txt[9:12]
                                                G = txt[12:15]
                                                X = txt[15:18]
                                            elif FS == 64:
                                                txt = txt064
                                                O = txt[9:12]
                                                G = txt[12:15]
                                                X = txt[15:18]
                                            elif FS == 128:
                                                txt = txt128
                                                O = txt[9:12]
                                                G = txt[12:15]
                                                X = txt[15:18]
                                            elif FS == 256:
                                                txt = txt256
                                                O = txt[9:12]
                                                G = txt[12:15]
                                                X = txt[15:18]

                                            O = np.reshape(O, (3, 1))

                                            S = np.array([[np.power((G[0] * 0.01), -1),
                                                           ((np.power((X[0] * 0.01 + 250), -1)) - 0.004),
                                                           (np.power((X[1] * 0.01 + 250), -1) - 0.004)],
                                                          [((np.power((X[0] * 0.01 + 250), -1)) - 0.004),
                                                           (np.power(G[1] * 0.01, -1)),
                                                           (np.power((X[2] * 0.01 + 250), -1) - 0.004)],
                                                          [(np.power((X[1] * 0.01 + 250), -1) - 0.004),
                                                           (np.power((X[2] * 0.01 + 250), -1) - 0.004),
                                                           (np.power(G[2] * 0.01, -1))]])

                                # Check for 0 length Record (i.e. USB connect) #
                                if 1 < len(record.payload) / 2:

                                # Check to see if date/time filter is enabled and within range #
                                #if NoFilter or (begin_timestamp < record.timestamp < end_timestamp):

                                    activity_count = activity_count + 1
                                    adxl.write('%s,%d,%s\n' % (timestamp, activityRecordLen, str(unixTime)))

                                    ##############################
                                    # Check for Timestamp issues #
                                    ##############################
                                    currenttimestamp = record.timestamp
                                    difference = currenttimestamp - lasttimestamp
                                    diff_in_s = difference.total_seconds()
                                    if (abs(diff_in_s) > 1 or diff_in_s < 0) and activity_count > 2:
                                        sumFile.write('%-35s %-15s %10s\n' % ('Timestamp Gap @ : ',
                                                                              str(record.timestamp.strftime(FMT)),
                                                                              str(unixTime)))
                                        timestampGap = timestampGap + 1
                                        fout1.write("%s,%s,%s\n" % (currenttimestamp, lasttimestamp, diff_in_s))
                                    lasttimestamp = currenttimestamp
                                    #################################################################################

                                    ##########################
                                    # Parse activity records #
                                    ##########################
                                    h = BitStream(hex=record.payload)
                                    lengthBytes = int(len(h.bytes))
                                    if record.type == 27 or record.type == 0:  # Nebula or Moses
                                        incrementsize += lengthBytes
                                        accelFS = lengthBytes / 4.5
                                    else:  # Taso
                                        incrementsize += lengthBytes
                                        accelFS = lengthBytes / 6
                                        isCalibrated = True

                                    if Log_Activity_Data:
                                        current_timedelta_accel = float(1 / accelFS)
                                        for j in range(0, int(accelFS)):
                                            timestamp = record.timestamp.strftime(FMT)
                                            if Log_Activity_Data:
                                                fout.write("%s," % timestamp)
                                                timedelta_accel += current_timedelta_accel
                                                fout.write("%.3f," % timedelta_accel)
                                            if record.type == 27:  # Nebula
                                                r = h.readlist("int:n, int:n, int:n", n=12)
                                                x, y, z = map(lambda f: f / 1, r)
                                                if isCalVals:
                                                    V = np.array([[x], [y], [z]])
                                                    accel = np.dot(S, (V - O))
                                                    x = float(accel[0][0])
                                                    y = float(accel[1][0])
                                                    z = float(accel[2][0])
                                            elif record.type == 0: # Moses
                                                r = h.readlist("int:n, int:n, int:n", n=12)
                                                x, y, z = map(lambda t: t / 256.0, r)
                                            else:  # Taso
                                                r = h.readlist("intle:n, intle:n, intle:n", n=16)
                                                x, y, z = map(lambda t: t / 256.0, r)

                                            numberRows = numberRows + 1
                                            vm = math.sqrt(x ** 2 + y ** 2 + z ** 2)
                                            fout.write("%.6f,%.6f,%.6f,%.6f\n" % (x, y, z, vm))
                                            ####################################################################################

                                            #########################
                                            # CHECK FOR PEGGED DATA #
                                            #########################
                                            if check_for_peg(x, y, z, activity_count, isCalibrated):
                                                if peggedState is False:
                                                    print('%-35s %-15s %10s' % ('peg @ : ', timestamp, str(unixTime)))
                                                    sumFile.write('%-35s %-15s %10s\n' % ('peg @ : ', timestamp,
                                                                                          str(unixTime)))
                                                    pegs = pegs + 1
                                                    peggedState = True
                                            else:
                                                if peggedState is True:
                                                    print('%-35s %-15s %10s' % (
                                                        'peg cleared @ : ', timestamp, str(unixTime)))
                                                    sumFile.write('%-35s %-15s %10s\n' % ('peg cleared @ : ', timestamp,
                                                                                          str(unixTime)))
                                                peggedState = False
                                            ####################################################################################

                                            ########################
                                            # Check For Flat Areas #
                                            ########################
                                            if check_for_flat(x, y, z, activity_count, isCalVals, accelFS):
                                                if flatRegionState is False and peggedState is False:
                                                    print('%-35s %-15s %10s' % ('Flat Region @ : ', timestamp,
                                                                                str(unixTime)))
                                                    sumFile.write('%-35s %-15s %10s\n' % ('Flat Region @ : ',
                                                                                          timestamp, str(unixTime)))
                                                    flatRegionState = True
                                                    flats = flats + 1
                                            else:
                                                if flatRegionState is True:
                                                    print('%-35s %-15s %10s' % ('Flat region cleared @ : ', timestamp,
                                                                                str(unixTime)))
                                                    sumFile.write(
                                                        '%-35s %-15s %10s\n' % ('Flat region cleared @ : ', timestamp,
                                                                                str(unixTime)))
                                                    flatRegionState = False
                                            ####################################################################################

                                else:
                                    timestamp = record.timestamp.strftime(FMT)
                                    print('%-35s %-15s %10s' % ('USB DOCK @ : ', timestamp, str(unixTime)))
                                    sumFile.write('%-35s %-15s %10s\n' % ('USB DOCK @ : ', timestamp, str(unixTime)))
                                i = i + 1

                            if LOG_IMU:
                                ###########################
                                # Parse STM32 IMU Schema  #
                                ###########################
                                #if record.type == 30 or record.type == 24:
                                if record.type == 24:
                                    # timestamp = record.timestamp.strftime(FMT)
                                    imuHeader = []
                                    isIMUdata = 1
                                    fout_imu = open(str(file) + '_imu.csv', "w")
                                    fout_imu.write("Timestamp,t,")
                                    h = BitStream(hex=record.payload)

                                    h.readlist("uintle:n", n=16)
                                    IMU_COLUMN = h.readlist("uintle:n", n=16)
                                    h.readlist("uintle:n", n=16)

                                    n = 0
                                    payload_length = 0
                                    while n < IMU_COLUMN[0]:
                                        IMU_COLUMN_FLAGS = h.readlist("uint:n", n=8)
                                        IMU_COLUMN_OFFSET = h.readlist("uint:n", n=8)
                                        IMU_COLUMN_SIZE = h.readlist("uint:n", n=8)
                                        IMU_SCALE_FACTOR = h.readlist("uint:n", n=32)
                                        IMU_COLUMN_HEADER = h.readlist("bytes:n", n=16)
                                        imuHeader.append(IMU_COLUMN_HEADER[0].decode('utf-8'))
                                        fout_imu.write("%s," % (IMU_COLUMN_HEADER[0].decode('utf-8')))
                                        print(IMU_COLUMN_FLAGS)
                                        print(IMU_COLUMN_OFFSET)
                                        print(IMU_COLUMN_SIZE)
                                        payload_length += (IMU_COLUMN_SIZE[0]/8)
                                        print(IMU_SCALE_FACTOR)
                                        print(IMU_COLUMN_HEADER[0].decode('utf-8'))
                                        n += 1
                                    fout_imu.write("\n")

                                # ###########################
                                # # Parse STM32 IMU Records #
                                # ###########################
                                # if record.type == 31:
                                #     #timestamp = record.timestamp.strftime(FMT)
                                #     h = BitStream(hex=record.payload)
                                #     FS2 = int(len(h.bytes) / payload_length)
                                #     current_timedelta = float(1/FS2)
                                #     for z in range(0, FS2):
                                #         #timestamp = record.timestamp.strftime(FMT)
                                #         timedelta += current_timedelta
                                #         fout_imu.write("%s," % timedelta)
                                #         imu = h.readlist("int:n, int:n, int:n, int:n, int:n, int:n", n=16)
                                #         temp = h.readlist("int:n", n=8)
                                #         ts = h.readlist("uint:n", n=16)
                                #         xa, ya, za, xg, yg, zg = map(lambda gy: gy / 1, imu)
                                #         temp1 = temp[0]
                                #         ts1 = ts[0]
                                #         fout_imu.write("%f,%f,%f,%f,%f,%f,%f,%f\n" % (xa, ya, za, xg, yg, zg, temp1, ts1))

                                ###########################
                                # Parse Taso IMU Records #
                                ###########################
                                if record.type == 25:
                                    # timestamp = record.timestamp.strftime(FMT)
                                    h = BitStream(hex=record.payload)
                                    FS2 = int(len(h.bytes) / payload_length)
                                    current_timedelta = float(1 / FS2)
                                    h.readlist("int:n", n=16)
                                    for p in range(0, FS2):
                                        timestamp = record.timestamp.strftime(FMT)
                                        timedelta += current_timedelta
                                        fout_imu.write("%s," % timestamp)
                                        fout_imu.write("%s," % timedelta)
                                        imu_acc = h.readlist("int:n, int:n, int:n", n=16)
                                        imu_temp = h.readlist("int:n", n=16)
                                        imu_gyr = h.readlist("int:n, int:n, int:n", n=16)
                                        imu_mag = h.readlist("int:n, int:n, int:n", n=16)
                                        h.readlist("int:n", n=8)
                                        xa, ya, za = map(lambda gy: float(gy) / 1, imu_acc)
                                        temperature = (float(imu_temp[0]) / 1) + 21
                                        xg, yg, zg = map(lambda gy: float(gy) / 1, imu_gyr)
                                        xm, ym, zm = map(lambda gy: float(gy) / 1, imu_mag)
                                        fout_imu.write(
                                            "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n" % (xa, ya, za, temperature, xg, yg, zg, xm, ym, zm))

                            ########################
                            # Parse Battery Record #
                            ########################
                            if record.type == 2:
                                batt = bytearray.fromhex(record.payload)
                                batt.reverse()
                                batt = ''.join(format(x, '02x') for x in batt)
                                timestamp = record.timestamp.strftime(FMT)
                                try:
                                    batt = int(batt, 16) * 0.001
                                except:
                                    print('issue with battery record')
                                fout2.write("%s,%f\n" % (timestamp, batt))
                                incrementsize += len(record.payload)/2

                            ########################
                            # Parse Temp Records #
                            ########################
                            if record.type == 30:
                                timedelta_temp += 4 # TODO Make adjustable
                                t = bytearray.fromhex(record.payload)
                                temp_sensor_type1 = t[:1]
                                stm32_temp = t[1:3]
                                stm32_temp.reverse()
                                stm32_temp = ''.join(format(x, '02x') for x in stm32_temp)
                                stm32_temp = struct.unpack('>H', bytes.fromhex(stm32_temp))

                                temp_sensor_type2 = t[3:4]
                                adxl_temp = t[4:6]
                                adxl_temp.reverse()
                                adxl_temp = ''.join(format(x, '02x') for x in adxl_temp)
                                adxl_temp = struct.unpack('>h', bytes.fromhex(adxl_temp))

                                timestamp = record.timestamp.strftime(FMT)

                                if temperatureCalibration:
                                    stm32_temp_cal = (gainMcu * stm32_temp[0]) + offsetMcu
                                    adxl_temp_cal = (gainAdxl * adxl_temp[0]) + offsetAdxl
                                    fout3.write("%s,%d,%d,%d,%d,%d,%0.2f,%0.2f\n" % (timestamp, timedelta_temp,
                                                                               temp_sensor_type1[0], stm32_temp[0],
                                                                               temp_sensor_type2[0], adxl_temp[0],
                                                                               stm32_temp_cal, adxl_temp_cal))
                                else:
                                    fout3.write("%s,%d,%d,%d,%d,%d\n" % (timestamp, timedelta_temp,
                                                                               temp_sensor_type1[0], stm32_temp[0],
                                                                               temp_sensor_type2[0], adxl_temp[0]))
                                incrementsize += len(record.payload) / 2

                            if record.type == 31:
                                timedelta_temp += 4  # TODO Make adjustable
                                t = bytearray.fromhex(record.payload)
                                temp_sensor_type1 = t[:1]
                                tmp117_temp = t[1:3]
                                tmp117_temp.reverse()
                                tmp117_temp = ''.join(format(x, '02x') for x in tmp117_temp)
                                tmp117_temp = struct.unpack('>H', bytes.fromhex(tmp117_temp))

                                timestamp = record.timestamp.strftime(FMT)

                                fout3.write("%s,%d,%d,%d\n" % (timestamp, timedelta_temp,
                                                                     temp_sensor_type1[0], tmp117_temp[0]))

                                incrementsize += len(record.payload) / 2


                            ###########################
                            # Parse Event Type Record #
                            ###########################
                            if record.type == 3:
                                incrementsize += (len(record.payload)/2)
                                timestamp = record.timestamp.strftime(FMT)
                                if record.payload == '0d':
                                    expected_reset_count = expected_reset_count + 1
                                    print('%-35s %-15s %10s' % ('Expected Reset @ :', timestamp, str(unixTime)))
                                    sumFile.write(
                                        '%-35s %-15s %10s\n' % ('Expected Reset @ :', timestamp, str(unixTime)))
                                elif record.payload == '01':
                                    unexpected_reset_count = unexpected_reset_count + 1
                                    print('%-35s %-15s %10s' % ('Unexpected Reset @ :', timestamp, str(unixTime)))
                                    sumFile.write(
                                        '%-35s %-15s %10s\n' % ('Unexpected Reset @ :', timestamp, str(unixTime)))
                                elif record.payload == '08':
                                    #print('%-35s %-15s %10s' % ('ENTER IDLE SLEEP @ ', timestamp, str(unixTime)))
                                    sumFile.write(
                                        '%-35s %-15s %10s\n' % ('ENTER IDLE SLEEP @ ', timestamp, str(unixTime)))
                                elif record.payload == '09':
                                    #print('%-35s %-15s %10s' % ('EXIT IDLE SLEEP @ ', timestamp, str(unixTime)))
                                    sumFile.write(
                                        '%-35s %-15s %10s\n' % ('EXIT IDLE SLEEP @ ', timestamp, str(unixTime)))
                                else:
                                    print('%-35s %-15s %10s' % (('Event Type %s @ ' % record.payload),
                                                                timestamp, str(unixTime)))
                                    sumFile.write('%-35s %-15s %10s\n' % (('Event Type %s @ ' % record.payload),
                                                                          timestamp, str(unixTime)))

                            ###########################
                            # Parse FIFO_ERROR Record #
                            ###########################
                            if record.type == 19:
                                incrementsize += (len(record.payload)/2)
                                timestamp = record.timestamp.strftime(FMT)
                                print('%-35s %-15s %10s' % ('FIFO ERROR @ : ', timestamp, str(unixTime)))
                                sumFile.write('%-35s %-15s %10s\n' % ('FIFO ERROR @ : ', timestamp, str(unixTime)))

                            ############################
                            # Parse DEBUG_ERROR Record #
                            ############################
                            if record.type == 23:
                                incrementsize += (len(record.payload)/2)
                                timestamp = record.timestamp.strftime(FMT)
                                print('%-35s %-15s %10s' % ('DEBUG_ERROR @ : ', timestamp, str(unixTime)))
                                sumFile.write('%-35s %-15s %10s\n' % ('DEBUG_ERROR @ : ', timestamp, str(unixTime)))

                            #########################
                            # Parse RAM_DUMP Record #
                            #########################
                            if record.type == 24:
                                incrementsize += (len(record.payload)/2)
                                timestamp = record.timestamp.strftime(FMT)
                                print('%-35s %-15s %10s' % ('RAM_DUMP @ : ', timestamp, str(unixTime)))
                                sumFile.write('%-35s %-15s %10s\n' % ('RAM_DUMP @ : ', timestamp, str(unixTime)))

                            #############################
                            # Parse EVENT_MARKER Record #
                            #############################
                            if record.type == 28:
                                incrementsize += (len(record.payload)/2)
                                timestamp = record.timestamp.strftime(FMT)
                                print('%-35s %-15s %10s' % ('EVENT MARKER @ : ', timestamp, str(unixTime)))
                                sumFile.write('%-35s %-15s %10s\n' % ('EVENT MARKER @ : ', timestamp, str(unixTime)))

                            ############################################
                            # Parse BUTTON_PRESS/BUTTON_RELEASE Record #
                            ############################################
                            if record.type == 29:
                                incrementsize += (len(record.payload)/2)
                                timestamp = record.timestamp.strftime(FMT)
                                if record.payload == '00':
                                    print('%-35s %-15s %10s' % ('BUTTON PRESS @ : ', timestamp, str(unixTime)))
                                    sumFile.write(
                                        '%-35s %-15s %10s\n' % ('BUTTON PRESS @ : ', timestamp, str(unixTime)))

                                if record.payload == '01':
                                    print('%-35s %-15s %10s' % ('BUTTON RELEASE @ : ', timestamp, str(unixTime)))
                                    sumFile.write(
                                        '%-35s %-15s %10s\n' % ('BUTTON RELEASE @ : ', timestamp, str(unixTime)))

                            ############################
                            # Parse ADXL_REGISTER_DUMP #
                            ############################
                            if record.type == 40:
                                incrementsize += (len(record.payload)/2)
                                timestamp = record.timestamp.strftime(FMT)
                                print(timestamp)
                                h = BitStream(hex=record.payload)
                                # s = Bits(h)
                                for j in range(0, 36):
                                    r = h.readlist("hex:n, hex:n", n=8)
                                    reg, value = map(lambda t: t, r)
                                    print(reg + ' ' + value)

                            ####################
                            # Parse Taso Epoch #
                            ####################
                            if record.type == 99: ###TODO
                                incrementsize += (len(record.payload) / 2)
                                timestamp = record.timestamp.strftime(FMT)
                                print(timestamp)
                                h = BitStream(hex=record.payload)
                                # s = Bits(h)
                                for j in range(0, 36):
                                    r = h.readlist("hex:n, hex:n", n=8)
                                    reg, value = map(lambda t: t, r)
                                    print(reg + ' ' + value)

                            ##########################
                            # Parse Calibration Data #
                            ##########################

                            if record.type in range(100, 114):
                                fout_cal.write("LOG_CAL_ORIENTATION_%s\n" % str(record.type - 100))
                                cal_orienatation = True

                            if record.type == 200 or record.type == 201 or record.type == 202 or record.type == 203:
                                if record.type == 200 and cal_orienatation is True:
                                    fout_cal.write("32HZ_CAL_DATA\n")
                                    cal_orienatation = False
                                if record.type == 201 and cal_orienatation is False:
                                    fout_cal.write("64HZ_CAL_DATA\n")
                                    cal_orienatation = True
                                if record.type == 202 and cal_orienatation is True:
                                    fout_cal.write("128HZ_CAL_DATA\n")
                                    cal_orienatation = False
                                if record.type == 203 and cal_orienatation is False:
                                    fout_cal.write("256HZ_CAL_DATA\n")
                                    cal_orienatation = True
                                h = BitStream(hex=record.payload)
                                #timestamp = record.timestamp.strftime(FMT)
                                #fout_cal.write("%s," % timestamp)
                                #timedelta_accel += current_timedelta_accel
                                #fout_cal.write("%.3f," % timedelta_accel)
                                cal = h.readlist("intle:n, intle:n, intle:n", n=16)
                                x, y, z = map(lambda t: t / 1.0, cal)
                                fout_cal.write("%.6f,%.6f,%.6f\n" % (x, y, z))

                            ##################################################
                            # Write Activity Data to file is option selected #
                            ##################################################
                            if Log_Activity_Data:
                                if activity_count > 1 and numberRows >= 1048320:
                                    numberRows = 0
                                    fout.close()
                                    if firstFile:
                                        unixTimeFileFirst = (firstTimestamp - datetime.fromtimestamp(0,timezone.utc)).total_seconds()
                                        try:
                                            os.rename(outPutPath + adxl_dir + filename + '.csv',
                                                      outPutPath + adxl_dir + filename + '-'
                                                      + str(int(unixTimeFileFirst)) + '.csv')
                                        except:
                                            os.remove(outPutPath + adxl_dir + filename + '-'
                                                      + str(int(unixTimeFileFirst)) + '.csv')
                                            os.rename(outPutPath + adxl_dir + filename + '.csv',
                                                      outPutPath + adxl_dir + filename + '-'
                                                      + str(int(unixTimeFileFirst)) + '.csv')
                                        firstFile = False
                                    unixTimeFileNext = (record.timestamp - datetime.fromtimestamp(0,
                                                                                                  timezone.utc)).total_seconds()
                                    fout = open(outPutPath + adxl_dir + filename
                                                + '-' + str(int(unixTimeFileNext)) + '.csv', "w")

                                    fout.write('ts,t,x,y,z,vm\n')


                        progress += (incrementsize + 9)
                        incrementsize = 0

                    unixTime = (lasttimestamp - datetime.fromtimestamp(0, timezone.utc)).total_seconds()
                    print('%-35s %-15s %10s' % ('Last Timestamp: ', lasttimestamp.strftime(FMT), str(unixTime)))
                    print('%-35s %-15s ' % ('Total Time: ', (str(unixTime - firstTimestampUnix))))
                    print(' ')
                    print('%-45s %5s' % (('Total number of Activity (%s) records: ' %
                                          hex(record_type)), str(activity_count)))
                    print('%-45s %5s' % ('Total number of Unexpected Resets: ', str(unexpected_reset_count)))
                    print('%-45s %5s' % ('Total number of Expected Resets: ', str(expected_reset_count)))
                    print('%-45s %5s' % ('Total number of Pegs: ', str(pegs)))
                    print('%-45s %5s' % ('Total flat areas: ', str(flats)))
                    print('%-45s %5s' % ('Timestamp Gaps: ', str(timestampGap)))
                    print(' ')
                    print('%-35s %-15s' % ('Parsing Start Date/Time: ', str(parse_start_date.strftime(FMT))))
                    print('%-35s %-15s' % ('Parsing Complete Date/Time: ', str(datetime.now().strftime(FMT))))
                    print(
                        '%-35s %-15s' % ('Total Time to Parse: ', str(datetime.now() - parse_start_date)))
                    print(' ')

                    sumFile.write('%-35s %-15s %10s\n\n' % ('Last Timestamp: ', lasttimestamp.strftime(FMT), str(unixTime)))
                    sumFile.write('%-45s %5s\n' % (('Total number of Activity (%s) records: ' %
                                                    hex(record_type)), str(activity_count)))
                    sumFile.write(
                        '%-45s %5s\n' % ('Total number of Unexpected Resets: ', str(unexpected_reset_count)))
                    sumFile.write('%-45s %5s\n' % ('Total number of Expected Resets: ', str(expected_reset_count)))
                    sumFile.write('%-45s %5s\n' % ('Total number of Pegs: ', str(pegs)))
                    sumFile.write('%-45s %5s\n' % ('Total flat areas: ', str(flats)))
                    sumFile.write('%-45s %5s\n' % ('Timestamp Gaps: ', str(timestampGap)))
                    sumFile.write('\n')
                    sumFile.write('%-35s %-15s\n' % ('Parsing Start Date/Time: ', str(parse_start_date.strftime(FMT))))
                    sumFile.write('%-35s %-15s\n' % ('Parsing Complete Date/Time: ', str(datetime.now().strftime(FMT))))
                    sumFile.write('%-35s %-15s\n' % ('Total Time to Parse: ',
                                                     str(datetime.now() - parse_start_date)))
                    sumFile.write('\n')
                    sumFile.flush()
                    os.fsync(sumFile)

                    fout_summary.write('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (str(file),firmware_version.rstrip(),
                                                                            str(firstTimestampUnix),str(unixTime),
                                                                            str(activity_count),str(unexpected_reset_count),
                                                                            str(expected_reset_count),str(pegs),str(flats),
                                                                            str(timestampGap)))

                    fout_summary.flush()
                    os.fsync(fout_summary)


                    if Log_Activity_Data:
                        fout.close()
                        # read in data
                        if createHtmlPlot:
                            df = pd.read_csv(str(file) + '.csv')
                            fig = px.line(df, x="t", y=['x', 'y', 'z', 'vm'], title='Acceleration',
                                          labels={"value": "acceleration in G",
                                                  "t": "Seconds"
                                                  }
                                          )
                            fig.write_html(str(file) + '.html')
                    if isIMUdata:
                        fout_imu.close()
                        df_imu = pd.read_csv(str(file) + '_imu.csv')
                        fig = px.line(df_imu, x="t", y=imuHeader[0:3], title='Acceleration',
                                      labels={"value": "acceleration in G",
                                              "Timestamp": "Seconds"
                                              }
                                      )
                        fig2 = px.line(df_imu, x="t", y=imuHeader[4:7], title='Gyro',
                                      labels={"value": "degrees/sec",
                                              "Timestamp": "Seconds"
                                              }
                                      )
                        fig3 = px.line(df_imu, x="t", y=imuHeader[3], title='Temperature',
                                       labels={imuHeader[6]: "Celsius",
                                               "Timestamp": "Seconds"
                                               }
                                       )
                        fig.write_html(str(file) + '_imu_accel.html')
                        fig2.write_html(str(file) + '_imu_gyro.html')
                        fig3.write_html(str(file) + '_imu_temp.html')
                    fout1.close()
                    fout2.close()
                    fout3.close()
                    fout_cal.close()
                    progress = 0
                if zipfiles:
                    try:
                        os.remove(basePath + "/log.bin")
                    except:
                        pass
                    try:
                        os.remove(basePath + "/calibration.json")
                    except:
                        pass
                    try:
                        os.remove(basePath + "/info.json")
                    except:
                        pass
                    try:
                        os.remove(basePath + "/temperature_calibration.json")
                    except:
                        pass
    sumFile.close()
    fout_summary.close()
    print('FINISHED\n')


def the_gui():

    global message, progress

    file_list_column = [
        [
            sg.Text("Image Folder"),
            sg.In(size=(25, 1), enable_events=True, key="-FOLDER-"),
            sg.FolderBrowse(),
        ],
        [sg.Listbox(values=[], enable_events=True, size=(40, 10), key="-FILE LIST-")],
        [
            sg.Button("PARSE FILE(S)", size=(20, 4)),
            sg.Checkbox("agdc/gt3x files", default=True, key="ZIPFILES"),
        ],
        [sg.Checkbox("Log Activity Data", default=False, key="LOGDATA")],
        [sg.Checkbox("Time/Date Filter", default=False, key="TIMEDATE")],
        [
            sg.Text("Beginning Unix Timestamp", size=(20, 1)),
            sg.InputText(key="first_timestamp", size=(15, 1)),
        ],
        [
            sg.Text("Ending Unix Timestamp", size=(20, 1)),
            sg.InputText(key="last_timestamp", size=(15, 1)),
        ],
        [sg.Checkbox("Apply Calibration if Possible", default=False, key="APPLYCAL")],
        [sg.Button("Exit")],
    ]

    image_viewer_column = [
        [sg.Output(size=(90, 30), background_color="black", text_color="white")]
    ]

    layout = [
        [
            sg.Column(file_list_column),
            sg.VSeperator(),
            sg.Column(image_viewer_column),
        ],
        [
            sg.Text('Work progress'),
            sg.ProgressBar(totallogsize, size=(20, 20), orientation='h', key='-PROG-')
        ]
    ]

    window = sg.Window("Nebula-Moses-Taso Parser: %s" % VERSION, layout)
    t = None

    while True:
        event, values = window.read(timeout=500)

        if event == sg.WIN_CLOSED or event == "Exit":
            break
            # Folder name was filled in, make a list of files in the folder
        if event == "-FOLDER-":
            folder = values["-FOLDER-"]
            try:
                # Get list of files in folder
                file_list = os.listdir(folder)
            except:
                file_list = []

            fnames = [
                f
                for f in file_list
                if os.path.isfile(os.path.join(folder, f)) and f.lower().endswith((".agdc", ".gt3x", ".bin", ".dat"))
            ]
            window["-FILE LIST-"].update(fnames)
        elif event == "-FILE LIST-":  # A file was chosen from the listbox
            try:
                filename = os.path.join(values["-FOLDER-"], values["-FILE LIST-"][0])
                window["-TOUT-"].update(filename)
                window["-IMAGE-"].update(filename=filename)
            except:
                pass
        if event == "PARSE FILE(S)":
            if values["LOGDATA"]:
                LOG_ACTIVITY_DATA = True
            else:
                LOG_ACTIVITY_DATA = False
            if values["TIMEDATE"]:
                NO_FILTER = 0
            else:
                NO_FILTER = 1
            if values["APPLYCAL"]:
                USE_CAL_VALUES = True
            else:
                USE_CAL_VALUES = False
            if values["ZIPFILES"]:
                USE_ZIP_FILES = True
            else:
                USE_ZIP_FILES = False

            t = Thread(target=main_process, args=(LOG_ACTIVITY_DATA, NO_FILTER, USE_CAL_VALUES, USE_ZIP_FILES,
                                                  values['-FOLDER-'],
                                                  values['first_timestamp'],
                                                  values['last_timestamp'], ), daemon=True)
            t.start()

        if t:
            window['-PROG-'].update_bar(progress, totallogsize)
            t.join(timeout=0)
            if not t.is_alive():                       # the thread finished
                # print(f'message = {message}')
                sg.popup_animated(None)                     # stop animination in case one is running
                t, message, progress = None, '', 0          # reset variables for next run
                window['-PROG-'].update_bar(0, 0)            # clear the progress bar


    window.close()


# Press the green button in the gutter to run the script.
if __name__ == "__main__":
    the_gui()
    print('Exiting Program')
