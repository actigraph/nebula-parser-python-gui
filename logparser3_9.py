from datetime import datetime, timezone
import collections, operator, struct
import functools


Record = collections.namedtuple('Record', 'type timestamp payload size bad_size')

def parse(fin, min_time, max_time):
    bad_size = 0
    try:
        while True:
            sync = ord(fin.read(1))
            while sync != 0x1E:
                sync = ord(fin.read(1))
                bad_size = bad_size + 1
            pos = fin.tell()
            buf = fin.read(7)
            chksum = functools.reduce(operator.xor, [(x) for x in buf], 0x1E)
            #chksum = 0x1E
            #for x in buf:
            #    chksum ^= x
            dtype, timestamp, size = struct.unpack('<BIH', buf)
            buf = fin.read(size)

            chksum = functools.reduce(operator.xor, [(x) for x in buf], chksum)
            checksum = (fin.read(1)).hex()
            if checksum == "":
                checksum = 0
            if chksum == 255 - int(checksum, 16):
                if (datetime.fromtimestamp(min_time, tz=timezone.utc)
                        < datetime.fromtimestamp(timestamp, tz=timezone.utc)
                        < datetime.fromtimestamp(max_time, tz=timezone.utc)):
                    yield Record(dtype, datetime.fromtimestamp(timestamp, tz=timezone.utc), buf.hex(), size, bad_size)
                else:
                    fin.seek(pos)
                    yield Record(253, datetime.fromtimestamp(timestamp, tz=timezone.utc), pos, size, bad_size)
                bad_size = 0
            else:
                # buf = pos
                fin.seek(pos)
                if (datetime.fromtimestamp(min_time, tz=timezone.utc)
                        < datetime.fromtimestamp(timestamp, tz=timezone.utc)
                        < datetime.fromtimestamp(max_time, tz=timezone.utc)):
                    yield Record(254, datetime.fromtimestamp(timestamp, tz=timezone.utc), pos, size, bad_size)
    except struct.error as e:
        return
    except TypeError as e:
        return


def datetime2timestamp(d):
    return int((d - datetime(1970, 1, 1)).total_seconds())


def pack(record):
    size = len(record.payload)
    fmt = '<BBIH{0}s'.format(size)
    t = datetime2timestamp(record.timestamp)
    buf = struct.pack(fmt, 0x1E, record.type, t, size, record.payload)
    chksum = 0
    for x in buf:
        chksum ^= ord(x)
    chksum = 255 - chksum
    buf += chr(chksum)
    return buf


def unpack(data):
    index = 0
    offset = 0
    sample = [0.0, 0.0, 0.0]

    while True:
        for i in range(3):
            if 0 == (offset & 0x7):
                if index == len(data):
                    return
                c = ord(data[index])
                index += 1
                shifter = ((c & 0xFF) << 4)
                offset += 8
                if index == len(data):
                    return
                c = ord(data[index])
                index += 1
                shifter |= ((c & 0xF0) >> 4)
                offset += 4
            else:
                shifter = ((c & 0x0F) << 8)
                offset += 4
                if index == len(data):
                    return
                c = ord(data[index])
                index += 1
                shifter |= (c & 0xFF)
                offset += 8

            # sign extern negative values
            if 0 != (shifter & 0x0800):
                shifter |= 0xF000 | -65536

            sample[i] = round(shifter / 256.0, 3)

        yield sample

