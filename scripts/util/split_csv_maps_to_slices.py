import argparse
import inspect
import logging
import os
import sys
from pathlib import Path

import coloredlogs
import numpy as np
import yaml

from pyspark import SparkContext, SparkConf, SQLContext
from pyspark.sql.types import *
from pyspark.ml.feature import Bucketizer
import pyspark.sql.functions as fn


_this_filename = inspect.getframeinfo(inspect.currentframe()).filename
_this_path = Path(_this_filename).parent.resolve()
sys.path.append(str(_this_path.parent))

import msproc

logging.basicConfig(format=msproc.LOG_FORMAT, level=logging.INFO)
logger = logging.getLogger(__name__)
coloredlogs.install(fmt=msproc.LOG_FORMAT, level='INFO', logger=logger)


def main():
    args = get_args()

    conf = SparkConf().setAppName(__name__)
    conf = conf.setMaster("local[*]")
    conf = conf.set('spark.local.dir', os.getenv('TMPDIR'))
    conf = conf.set('spark.executor.memory', '2500M').\
        set('spark.driver.memory', '10G')

    sc = SparkContext(conf=conf)
    sqlc = SQLContext(sc)
    logger.info('Got a Spark context')

    input_dir = Path(args.config['root_dir']) / args.config['samples_adjusted_swaths']
    input_filenames = [str(fname) for fname in input_dir.glob('*.csv')]

    if not input_filenames:
        raise Exception('No input files found at configured location: %s' % str(input_dir))
    else:
        logger.info('List of maps to split: \n %s' % ('\n'.join(input_filenames)))

    output_dir = str(Path(args.config['root_dir']) / args.config['slices_location'])

    ms_schema = StructType([StructField("spectrum_index", IntegerType()),
                            StructField("level", ByteType()),
                            StructField("rt", FloatType()),
                            StructField("mz", FloatType()),
                            StructField("intensity", FloatType()),
                            StructField("prec_mz", FloatType()),
                            StructField("swath_upper_adjusted", FloatType()),
                            StructField("swath_lower_adjusted", FloatType())])
    sample = sqlc.read.csv(input_filenames, schema=ms_schema, header="true")

    max_rt = sample.select('rt').rdd.max()[0]
    save_max_rt_to_computed_values_file(args, max_rt)
    logger.info('Read sample files')

    # ## Mark each time value with its corresponding window number

    # Construct Bucketizer for binning rt values
    rt_windows = np.arange(0, max_rt + args.config['window_size_sec'],
                           args.config['window_size_sec'])
    bucketizer = Bucketizer(splits=rt_windows, inputCol="rt", outputCol="rt_window")

    # Read, bin, and concatenate all input files
    samples = []
    for filename in input_filenames:
        sample = sqlc.read.csv(filename, schema=ms_schema, header="true")
        # Bin rt values
        sample = bucketizer.setHandleInvalid("keep").transform(sample)
        name = os.path.splitext(os.path.basename(filename))[0]
        sample = sample.withColumn("file", fn.lit(name))
        sample = sample.withColumn("swath_lower_adjusted",
                                   fn.regexp_replace(fn.format_number(sample['swath_lower_adjusted'], 2), ",", ""))
        sample = sample.withColumn("swath_upper_adjusted",
                                   fn.regexp_replace(fn.format_number(sample['swath_upper_adjusted'], 2), ",", ""))
        samples.append(sample.rdd)
    samples = sc.union(samples)
    logger.info('Binned sample values in RT windows')

    # ## Split according to swath and time window, and save to CSVs
    # Note:
    # May fail with `coalesce()` on too much data. Check memory available to driver.
    # It also makes the task slower.
    logger.info('Starting Spark execution and writing to CSVs ...')
    samples.toDF(sampleRatio=0.2) \
        .repartition(fn.col("swath_lower_adjusted"), fn.col("rt_window"), fn.col("file")) \
        .coalesce(1) \
        .write.format('com.databricks.spark.csv')\
        .mode('overwrite')\
        .partitionBy('swath_lower_adjusted', 'rt_window')\
        .option('header', 'true').csv(output_dir)
    logger.info('... done.')


def save_max_rt_to_computed_values_file(args, max_rt):
    values_filename =  Path(args.config['root_dir']) / 'computed_values.yaml'

    if not values_filename.exists():
        computed_values = {'max_rt': max_rt}
    else:
        with values_filename.open(mode='r') as values_file:
            computed_values = yaml.load(values_file)
        computed_values['max_rt'] = max_rt

    with values_filename.open(mode='w') as values_file:
        yaml.dump(computed_values, values_file)


def get_args():
    desc = 'Split SWATH maps as CSVs into (time window, swath) slices'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-c',
                        '--config',
                        required=True,
                        type=str,
                        help='YAML experiment config file')
    args = parser.parse_args()
    with open(args.config, 'r') as config_file:
        args.config = yaml.load(config_file)
    return args


if __name__ == '__main__':
    main()
