from contextlib import contextmanager, redirect_stderr, redirect_stdout
from os import devnull

import subprocess
from subprocess import DEVNULL

@contextmanager
def suppress_output(verbose):
    """Suppress output when """
    if verbose:
        pass
    else:
        with open(devnull, 'w') as fnull:
            with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
                yield (err, out)

def run_command(command, verbose):
    if verbose:
        subprocess.run(command, shell=True)
    else:
        subprocess.run(command, shell=True, stdout=DEVNULL, stderr=DEVNULL)

