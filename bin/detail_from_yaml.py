#!/usr/bin/env python3
# Daniel Stribling
# ORCID: 0000-0002-0649-9506
# Renne Lab, University of Florida
# hkref (hybkit Ref)  http://www.github.com/RenneLab/hkref
# Hybkit Project : http://www.github.com/RenneLab/hybkit

"""
Get a detail from a yaml file
"""

import argparse
import yaml

def detail_from_yaml(yaml_file, detail):
    # ---- Read Yaml File ----
    with open(yaml_file, "r") as yaml_file_obj:
        yaml_data = yaml.safe_load(yaml_file_obj)

    use_obj = yaml_data
    for key in detail:
        use_obj = use_obj[key]
    print(use_obj)


def parse_args(print_help=False):
    parser = argparse.ArgumentParser(description=__doc__)
    #                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--yaml_file', action='store', required=True,
                        help='YAML file to read')
    parser.add_argument('--detail', action='store', required=True,
                        nargs='+',
                        help='Space-separated list of mapping keys of detail to retreive')

    if print_help:
        parser.print_help()
        return {}

    namespace = parser.parse_args()
    args_info = vars(namespace)
    if not args_info:
        sys.exit()

    return args_info

if __name__ == '__main__':
    args_info = parse_args()
    detail_from_yaml(**args_info)

