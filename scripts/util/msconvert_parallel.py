# /usr/bin/python
import argparse
import os
import time
from multiprocessing import Pool
from subprocess import Popen, DEVNULL


parser = argparse.ArgumentParser(description='Runs msconvert in parallel on Windows')
parser.add_argument('-i', '--input_dir', required=True, type=str, help='Input directory')
parser.add_argument('-o', '--output_dir', default='./converted', type=str,
                    help='Ending added to input_dir for output')
parser.add_argument('-e', '--extension', default="wiff", type=str,
                    help="Input extension")
parser.add_argument('-n', '--number_threads', default=8, type=int,
                    help='Number of threads to use at a time')
parser.add_argument('-f', '--format', default="mzML", type=str, help='Format to convert to')
args = parser.parse_args()

input_dir = args.input_dir
output_dir = args.output_dir
extension = args.extension
new_format = args.format
nthreads = args.number_threads

msconvert_bin = "msconvert.exe"
msconvert_params = '--32 --zlib --noindex --filter "peakPicking vendor msLevel=1-2"'


def run_msconvert(input_filename):
    output_filename = os.path.splitext(input_filename)[0] + '.' + new_format
    command = " ".join([msconvert_bin,
                        "--" + new_format,
                        msconvert_params,
                        '-o', output_dir,
                        os.path.join(input_dir, input_filename),
                        "--outfile", os.path.join(output_dir, output_filename)])
    print(command)
    Popen(command, stdout=DEVNULL).wait()
    print("Wrote " + os.path.join(output_dir, output_filename))


if __name__ == '__main__':
    """
    Main process needs to be separated in own function
    to avoid recursive calls
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    input_filenames = list(filter(lambda fname: fname.endswith(extension),
                                  os.listdir(input_dir)))
    start_time = time.time()

    with Pool(processes=nthreads) as pool:
        pool.map(run_msconvert, input_filenames)

    elapsed_time = time.time() - start_time
    print("Total time {}".format(elapsed_time))
    os.system('pause')
