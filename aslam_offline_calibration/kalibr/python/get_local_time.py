"""
Load data with local and remote time, get the corrected local time synced to the remote time, and save the timestamps.
"""
import os

import numpy as np

import sm
import argparse


def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('inlog',
        help= "data log file in csv format. Each row has local host time and remote device time in seconds")
    parser.add_argument(
        '--localtimeindex',
        type=int, default=1,
        help="1-based column index of the local host time in the log. (default: %(default)s)"
        )
    parser.add_argument(
        '--remotetimeindex',
        type=int, default=2,
        help="1-based column index of the local host time in the log. (default: %(default)s)"
        )
    parser.add_argument(
        '--outlog',
        help="output log file. (default: %(default)s)")
    return parser.parse_args()


def main():
    args = parseArgs()
    data = np.loadtxt(args.inlog, delimiter=',')
    remotetimes = data[:, args.remotetimeindex - 1]
    localtimes = data[:, args.localtimeindex - 1]

    outputlog = args.outlog
    if not args.outlog:
        inlognoext = os.path.splitext(args.inlog)[0]
        outputlog = inlognoext + "-syncedlocaltimes.log"

    timestamp_corrector = sm.DoubleTimestampCorrector()
    for i, remotetime in enumerate(remotetimes):
        timestamp_corrector.correctTimestamp(remotetime, localtimes[i])
    correctedtimes = []
    for i, remotetime in enumerate(remotetimes):
        correctedtimes.append(timestamp_corrector.getLocalTime(remotetime))

    np.savetxt(outputlog, correctedtimes, fmt="%.9f", delimiter=",")
    print('Saved corrected local time in {}'.format(outputlog))


if __name__ == '__main__':
    main()
