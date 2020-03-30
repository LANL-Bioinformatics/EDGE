import psutil
import argparse
import os
import logging
from send_notifications import send_sms_message, send_email_message

import logging
logging.basicConfig(format='%(asctime)s %(message)s', filename='disk_usage.log',level=logging.DEBUG)

def bytes2human(n):
    # http://code.activestate.com/recipes/578019
    # >>> bytes2human(10000)
    # '9.8K'
    # >>> bytes2human(100001221)
    # '95.4M'
    symbols = ('K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y')
    prefix = {}
    for i, s in enumerate(symbols):
        prefix[s] = 1 << (i + 1) * 10
    for s in reversed(symbols):
        if n >= prefix[s]:
            value = float(n) / prefix[s]
            return '%.1f%s' % (value, s)
    return "%sB" % n

def disk_usage(vol, threshold):
    total, used, free, percent = psutil.disk_usage(vol)
    msg = 'total: {0} used: {1}'.format(bytes2human(total), bytes2human(used))
    logging.info(msg)
    if percent > threshold:
        msg = 'Percent Disk usage: {0}, remaining disk space: {1}'.format(str(percent), bytes2human(free))
        logging.info(msg)
        send_sms_message(msg)
        send_email_message(msg)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--volume', required=False, help='Disk volume to monitor')
    parser.add_argument('-t', '--threshold', required=False, help='Threshold for notification')
    args = vars((parser.parse_args()))
    volume = args['volume'] if args['volume'] else os.getcwd()
    threshold = float(args['threshold']) if args['threshold'] else 0.85
    logging.info('disk usage running on volume {0} with alert notification threshold {1}'.format(volume, threshold))
    disk_usage(vol=volume, threshold=threshold)