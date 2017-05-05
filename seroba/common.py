import os
import time
import sys
import subprocess
import urllib.request
import pyfastaq


class Error (Exception): pass


def syscall(cmd, allow_fail=False, verbose=False, verbose_filehandle=sys.stdout, print_errors=True):
    if verbose:
        print('syscall:', cmd, flush=True, file=verbose_filehandle)
    try:
        subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as error:
        errors = error.output.decode()
        if print_errors:
            print('The following command failed with exit code', error.returncode, file=sys.stderr)
            print(cmd, file=sys.stderr)
            print('\nThe output was:\n', file=sys.stderr)
            print(errors, file=sys.stderr, flush=True)

        if allow_fail:
            return False, errors
        else:
            sys.exit(1)

    return True, None


def decode(x):
    try:
        s = x.decode()
    except:
        return x
    return s

def detect_sequence_format(infile):

    sequence_file = pyfastaq.utils.open_file_read(infile)
    first = sequence_file.readline()
    pyfastaq.utils.close(sequence_file)
    if first[0] == '>':
        return 'fa'
    elif first[0] == '@':
        return 'fq'
    else:
        raise Error('wrong format of sequence_file. First character is '+first[0] )
