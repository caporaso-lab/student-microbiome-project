#!/usr/bin/env python
# File created on 22 Apr 2013
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from numpy import median
from cogent.maths.stats.test import correlation_test
from qiime.util import parse_command_line_parameters, make_option
from qiime.parse import parse_distmat, parse_mapping_file, parse_mapping_file_to_dict
from qiime.group import get_grouped_distances
from qiime.util import qiime_open
from qiime.filter import filter_samples_from_distance_matrix, sample_ids_from_metadata_description

script_info = {}
script_info['help_on_no_arguments'] = False
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = [("","","")]
script_info['output_description']= ""
script_info['required_options'] = []
script_info['optional_options'] = []
script_info['version'] = __version__

wdm_fp = "/Users/caporaso/analysis/student-microbiome-project/beta-diversity/weighted_unifrac_dm.txt.gz"
udm_fp = "/Users/caporaso/analysis/student-microbiome-project/beta-diversity/unweighted_unifrac_dm.txt.gz"
mapping_fp = './SMP-map.tsv'

def correlated_variability(grouped_distances1,
                           grouped_distances2,
                           summarize_f=median):
    pids1 = dict([(e[0],e[2]) for e in grouped_distances1])
    pids2 = dict([(e[0],e[2]) for e in grouped_distances2])
    
    common_pids = list(set(pids1.keys()) & set(pids2.keys()))
    d1 = [summarize_f(pids1[pid]) for pid in common_pids]
    d2 = [summarize_f(pids2[pid]) for pid in common_pids]
    r = correlation_test(d1,d2,method='spearman')
    return len(common_pids), r

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    md, mh, _ = parse_mapping_file(open(mapping_fp))
    body_sites = ['Gut','Tongue','Palm','Forehead']
    intraindividual_distances = []
    
    print "Unweighted UniFrac"
    for b in body_sites:
        dm_fp = "/Users/caporaso/analysis/student-microbiome-project/beta-diversity/unweighted_unifrac_dm.%s_ts_only.txt.gz" % b.lower()
        h, d = parse_distmat(qiime_open(dm_fp))
        intraindividual_distances.append(get_grouped_distances(h, d, mh, md, 'PersonalID'))
    
    for i in range(len(body_sites)):
        for j in range(i):
            r = correlated_variability(intraindividual_distances[i],intraindividual_distances[j])
            print "%s/%s (n=%d): rho:%1.3f, p=%f" % (body_sites[i],body_sites[j],r[0],r[1][0],r[1][3])
    
    print "**"
    print "Weighted UniFrac"
    for b in body_sites:
        dm_fp = "/Users/caporaso/analysis/student-microbiome-project/beta-diversity/weighted_unifrac_dm.%s_ts_only.txt.gz" % b.lower()
        h, d = parse_distmat(qiime_open(dm_fp))
        intraindividual_distances.append(get_grouped_distances(h, d, mh, md, 'PersonalID'))
    
    for i in range(len(body_sites)):
        for j in range(i):
            r = correlated_variability(intraindividual_distances[i],intraindividual_distances[j])
            print "%s/%s (n=%d): rho:%1.3f, p=%f" % (body_sites[i],body_sites[j],r[0],r[1][0],r[1][3])



if __name__ == "__main__":
    main()