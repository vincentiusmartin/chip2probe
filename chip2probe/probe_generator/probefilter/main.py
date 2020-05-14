"""
This is the main file used in running probefilter.

The output is a file of clean non custom probes as well as the optional
clean custom probes, with primers

Authors: Farica Zhuang, Vincentius Martin
Created on Sep 30, 2019
"""

import yaml
import sys

#read the config file here
with open("../config.yml", "r") as ymlfile:
    conf = yaml.load(ymlfile, Loader=yaml.FullLoader)['default']

from probefilter import ProbeFilter

if __name__ == '__main__':
    if len(sys.argv) < 2:
        raise Exception("Please input sequence list")
    input_data = sys.argv[1]

    ProbeFilter(input_data, **conf)._run_all()
