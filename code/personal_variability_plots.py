#!/usr/bin/env python
# File created on 15 Mar 2013
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"


from qiime.util import parse_command_line_parameters, make_option
from pylab import hist, xlim, figure, plot, subplot, savefig, ylim
from qiime.group import get_adjacent_distances
from cogent.draw.distribution_plots import generate_box_plots
from qiime.parse import parse_distmat, parse_mapping_file_to_dict
from qiime.util import qiime_open
from cogent.maths.stats.test import mc_t_two_sample

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = [("","","")]
script_info['output_description']= ""
script_info['required_options'] = [\
 # Example required option
 make_option('-m','--mapping_fp',type="existing_filepath",help='the input mapping file'),\
]
script_info['optional_options'] = [
 make_option('-o','--output_dir',type="new_dirpath",
             help='the output directory [default: %default]',default='./'),\
]
script_info['version'] = __version__

wdm_fp = "/Users/caporaso/analysis/student-microbiome-project/beta-diversity/weighted_unifrac_dm.txt.gz"
udm_fp = "/Users/caporaso/analysis/student-microbiome-project/beta-diversity/unweighted_unifrac_dm.txt.gz"



def rank_test(data,value='Yes',tails='high'):
    data.sort()
    v_dist = []
    non_v_dist = []
    for d in data:
        if d[1] == value:
            v_dist.append(d[0])
        else:
            non_v_dist.append(d[0])
    t, parametric_p, ts, non_parametric_p = mc_t_two_sample(v_dist,non_v_dist,tails=tails)
    return t, parametric_p, non_parametric_p, v_dist, non_v_dist

def rank_adjacent_unifrac(sid_data,
                          dm_header,
                          dm_data,
                          map_data,
                          time_field = 'WeeksSinceStart',
                          line_color='blue'):
    disturbed_weeks = [float(e[0]) for e in sid_data if e[2]]
    if len(disturbed_weeks) < 1:
        raise ValueError, "Must be some disturbance events."
    adjacent_distances = get_adjacent_distances(dm_header,dm_data,[e[1] for e in sid_data])
    distance_to_disturbance = []
    for d, sample_ids in zip(adjacent_distances[0],adjacent_distances[1]):
        distance_to_disturbance.append((d,map_data[sample_ids[1]]['SampleAntibioticDisturbance']))
    
    distance_to_disturbance.sort()
    result = []
    for i, dd in enumerate(distance_to_disturbance):
        result.append((i,dd[1],dd[0]))
    return result

def score_ranked_adjacent_unifracs(mapping_data,
                                   dm,
                                   dm_header,
                                   inclusion_field,
                                   inclusion_value,
                                   personal_id_field,
                                   disturbed_field,
                                   disturbed_value,
                                   time_field="WeeksSinceStart"):
    sids = get_timeseries_sample_ids(mapping_data,
                                     inclusion_field=inclusion_field,
                                     inclusion_value=inclusion_value,
                                     personal_id_field=personal_id_field,
                                     disturbed_fields=[disturbed_field],
                                     disturbed_value=disturbed_value,
                                     time_field=time_field)
    ranked_data = []
    for sid_data in sids.values():
        try:
            x = rank_adjacent_unifrac(sid_data,dm_header,dm,mapping_data)
        except ValueError:
            pass
        else:
            ranked_data.extend(x)

    test_result = rank_test(ranked_data)
    return test_result[1], test_result[2], test_result[3], test_result[4]

def get_timeseries_sample_ids(mapping_d,
                        inclusion_field,
                        inclusion_value,
                        personal_id_field,
                        disturbed_fields,
                        disturbed_value,
                        time_field):
    pid_to_adjacent_unifrac = {}
    for k, v in mapping_d.items():
        if v[inclusion_field] == inclusion_value:
            personal_id = v[personal_id_field]
            week = v[time_field]
            sample_disturbances = []
            for disturbed_field in disturbed_fields:
                if v[disturbed_field] == disturbed_value:
                    sample_disturbances.append(disturbed_field)
            try:
                pid_to_adjacent_unifrac[personal_id].append((week,k,sample_disturbances))
            except KeyError:
                pid_to_adjacent_unifrac[personal_id] = [(week,k,sample_disturbances)]
    
    for personal_id in pid_to_adjacent_unifrac:
        pid_to_adjacent_unifrac[personal_id].sort()
    
    return pid_to_adjacent_unifrac

def plot_adjacent_unifrac(sid_data,
                          dm_header,
                          dm_data,
                          map_data,
                          time_field = 'WeeksSinceStart',
                          line_color='blue',
                          disturbance_colors=['blue','green','red','orange']):
    disturbed_weeks = []
    disturbance_types = []
    for e in sid_data:
        for d in e[2]:
            disturbed_weeks.append((float(e[0]),d))
            disturbance_types.append(d)
    disturbance_types = list(set(disturbance_types))
    disturbance_colors = dict(zip(disturbance_types,disturbance_colors))
    adjacent_distances = get_adjacent_distances(dm_header,dm_data,[e[1] for e in sid_data])
    
    x = map(float,[map_data[e[1]][time_field] for e in adjacent_distances[1]])
    y = adjacent_distances[0]
    plot(x,y,'-.',c=line_color)
    ylim(0,1)
    ymin, ymax = ylim()
    for dw,dt in disturbed_weeks:
        line_color = disturbance_colors[dt]
        dw_w_offset = dw + (disturbance_types.index(dt) * 0.1)
        plot([dw_w_offset,dw_w_offset],[ymin,ymax],'--',c=line_color)

def plot_adjacent_unifracs(dm_header,
                           dm_data,
                           map_data,
                           inclusion_fields,
                           inclusion_value,
                           disturbed_fields,
                           disturbed_value,
                           output_fp,
                           personal_id_field = 'PersonalID',
                           time_field = 'WeeksSinceStart',
                           colors=['brown','red','blue','orange']):
    data = {}
    all_pids = []
    for inclusion_field in inclusion_fields:
        pid_data = get_timeseries_sample_ids(
                                         map_data,
                                         inclusion_field,
                                         inclusion_value,
                                         personal_id_field,
                                         disturbed_fields,
                                         disturbed_value,
                                         time_field)
        all_pids.extend(pid_data)
        data[inclusion_field] = pid_data
    all_pids = set(all_pids)
    
    # configure the plot and subplots
    num_rows = len(all_pids)
    num_cols = 4
    figure(num=None, figsize=(20, 200), dpi=80, facecolor='w', edgecolor='k')

    for i,pid in enumerate(all_pids):
        plot_num = i+1
        subplot(num_rows,num_cols,plot_num)
        for inclusion_field,color in zip(inclusion_fields,colors):
            try:
                datum = data[inclusion_field][pid]
            except KeyError:
                pass
            else:
                plot_adjacent_unifrac(datum,dm_header,dm_data,map_data,line_color=color)
    savefig(output_fp)

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)

    m = parse_mapping_file_to_dict(open(opts.mapping_fp,'U'))[0]
    wh, wdm = parse_distmat(qiime_open(wdm_fp,'U'))
    uh, udm = parse_distmat(qiime_open(udm_fp,'U'))

    plot_adjacent_unifracs(uh,udm,m,["GutTimeseries","TongueTimeseries","ForeheadTimeseries","PalmTimeseries"],
                           "Yes",["SampleAntibioticDisturbance","SampleMenstruationDisturbance","SampleSicknessDisturbance"],
                           "Yes",output_fp='unweighted-abx.pdf')
    plot_adjacent_unifracs(wh,wdm,m,["GutTimeseries","TongueTimeseries","ForeheadTimeseries","PalmTimeseries"],
                           "Yes",["SampleAntibioticDisturbance","SampleMenstruationDisturbance","SampleSicknessDisturbance"],
                           "Yes",output_fp='weighted-abx.pdf')

    r = score_ranked_adjacent_unifracs(m,
                                       udm,
                                       uh,
                                       inclusion_field="GutTimeseries",
                                       inclusion_value="Yes",
                                       personal_id_field="PersonalID",
                                       disturbed_field="SampleAntibioticDisturbance",
                                       disturbed_value="Yes")
    print "Unweighted UniFrac: %1.3f (n-disturbed samples: %d, n-undisturbed samples: %d) " % (r[1], len(r[2]), len(r[3])) 
    r = score_ranked_adjacent_unifracs(m,
                                       wdm,
                                       wh,
                                       inclusion_field="GutTimeseries",
                                       inclusion_value="Yes",
                                       personal_id_field="PersonalID",
                                       disturbed_field="SampleAntibioticDisturbance",
                                       disturbed_value="Yes")
    print "Weighted UniFrac: %1.3f (n-disturbed samples: %d, n-undisturbed samples: %d) " % (r[1], len(r[2]), len(r[3]))

    r = score_ranked_adjacent_unifracs(m,
                                       udm,
                                       uh,
                                       inclusion_field="TongueTimeseries",
                                       inclusion_value="Yes",
                                       personal_id_field="PersonalID",
                                       disturbed_field="SampleAntibioticDisturbance",
                                       disturbed_value="Yes")
    print "Unweighted UniFrac: %1.3f (n-disturbed samples: %d, n-undisturbed samples: %d) " % (r[1], len(r[2]), len(r[3]))
    r = score_ranked_adjacent_unifracs(m,
                                       wdm,
                                       wh,
                                       inclusion_field="TongueTimeseries",
                                       inclusion_value="Yes",
                                       personal_id_field="PersonalID",
                                       disturbed_field="SampleAntibioticDisturbance",
                                       disturbed_value="Yes")
    print "Weighted UniFrac: %1.3f (n-disturbed samples: %d, n-undisturbed samples: %d) " % (r[1], len(r[2]), len(r[3]))

    r = score_ranked_adjacent_unifracs(m,
                                       udm,
                                       uh,
                                       inclusion_field="PalmTimeseries",
                                       inclusion_value="Yes",
                                       personal_id_field="PersonalID",
                                       disturbed_field="SampleAntibioticDisturbance",
                                       disturbed_value="Yes")
    print "Unweighted UniFrac: %1.3f (n-disturbed samples: %d, n-undisturbed samples: %d) " % (r[1], len(r[2]), len(r[3]))
    r = score_ranked_adjacent_unifracs(m,
                                       wdm,
                                       wh,
                                       inclusion_field="PalmTimeseries",
                                       inclusion_value="Yes",
                                       personal_id_field="PersonalID",
                                       disturbed_field="SampleAntibioticDisturbance",
                                       disturbed_value="Yes")
    print "Weighted UniFrac: %1.3f (n-disturbed samples: %d, n-undisturbed samples: %d) " % (r[1], len(r[2]), len(r[3]))



    r = score_ranked_adjacent_unifracs(m,
                                       udm,
                                       uh,
                                       inclusion_field="ForeheadTimeseries",
                                       inclusion_value="Yes",
                                       personal_id_field="PersonalID",
                                       disturbed_field="SampleAntibioticDisturbance",
                                       disturbed_value="Yes")
    print "Unweighted UniFrac: %1.3f (n-disturbed samples: %d, n-undisturbed samples: %d) " % (r[1], len(r[2]), len(r[3]))
    r = score_ranked_adjacent_unifracs(m,
                                       wdm,
                                       wh,
                                       inclusion_field="ForeheadTimeseries",
                                       inclusion_value="Yes",
                                       personal_id_field="PersonalID",
                                       disturbed_field="SampleAntibioticDisturbance",
                                       disturbed_value="Yes")
    print "Weighted UniFrac: %1.3f (n-disturbed samples: %d, n-undisturbed samples: %d) " % (r[1], len(r[2]), len(r[3]))

    # 
    # r = score_ranked_adjacent_unifracs(m,
    #                                    udm,
    #                                    uh,
    #                                    inclusion_field="GutTimeseries",
    #                                    inclusion_value="Yes",
    #                                    personal_id_field="PersonalID",
    #                                    disturbed_field="SampleMenstruationDisturbance",
    #                                    disturbed_value="Yes")
    # print "Unweighted UniFrac: %1.3f (n-disturbed samples: %d, n-undisturbed samples: %d) " % (r[1], len(r[2]), len(r[3]))
    # r = score_ranked_adjacent_unifracs(m,
    #                                    wdm,
    #                                    wh,
    #                                    inclusion_field="GutTimeseries",
    #                                    inclusion_value="Yes",
    #                                    personal_id_field="PersonalID",
    #                                    disturbed_field="SampleMenstruationDisturbance",
    #                                    disturbed_value="Yes")
    # print "Weighted UniFrac: %1.3f (n-disturbed samples: %d, n-undisturbed samples: %d) " % (r[1], len(r[2]), len(r[3]))
    # 
    # # <codecell>
    # 
    # r = score_ranked_adjacent_unifracs(m,
    #                                    udm,
    #                                    uh,
    #                                    inclusion_field="TongueTimeseries",
    #                                    inclusion_value="Yes",
    #                                    personal_id_field="PersonalID",
    #                                    disturbed_field="SampleMenstruationDisturbance",
    #                                    disturbed_value="Yes")
    # print "Unweighted UniFrac: %1.3f (n-disturbed samples: %d, n-undisturbed samples: %d) " % (r[1], len(r[2]), len(r[3]))
    # r = score_ranked_adjacent_unifracs(m,
    #                                    wdm,
    #                                    wh,
    #                                    inclusion_field="TongueTimeseries",
    #                                    inclusion_value="Yes",
    #                                    personal_id_field="PersonalID",
    #                                    disturbed_field="SampleMenstruationDisturbance",
    #                                    disturbed_value="Yes")
    # print "Weighted UniFrac: %1.3f (n-disturbed samples: %d, n-undisturbed samples: %d) " % (r[1], len(r[2]), len(r[3]))
    # 
    # # <codecell>
    # 
    # r = score_ranked_adjacent_unifracs(m,
    #                                    udm,
    #                                    uh,
    #                                    inclusion_field="PalmTimeseries",
    #                                    inclusion_value="Yes",
    #                                    personal_id_field="PersonalID",
    #                                    disturbed_field="SampleMenstruationDisturbance",
    #                                    disturbed_value="Yes")
    # print "Unweighted UniFrac: %1.3f (n-disturbed samples: %d, n-undisturbed samples: %d) " % (r[1], len(r[2]), len(r[3]))
    # r = score_ranked_adjacent_unifracs(m,
    #                                    wdm,
    #                                    wh,
    #                                    inclusion_field="PalmTimeseries",
    #                                    inclusion_value="Yes",
    #                                    personal_id_field="PersonalID",
    #                                    disturbed_field="SampleMenstruationDisturbance",
    #                                    disturbed_value="Yes")
    # print "Weighted UniFrac: %1.3f (n-disturbed samples: %d, n-undisturbed samples: %d) " % (r[1], len(r[2]), len(r[3]))
    # 
    # # <codecell>
    # 
    # r = score_ranked_adjacent_unifracs(m,
    #                                    udm,
    #                                    uh,
    #                                    inclusion_field="ForeheadTimeseries",
    #                                    inclusion_value="Yes",
    #                                    personal_id_field="PersonalID",
    #                                    disturbed_field="SampleMenstruationDisturbance",
    #                                    disturbed_value="Yes")
    # print "Unweighted UniFrac: %1.3f (n-disturbed samples: %d, n-undisturbed samples: %d) " % (r[1], len(r[2]), len(r[3]))
    # r = score_ranked_adjacent_unifracs(m,
    #                                    wdm,
    #                                    wh,
    #                                    inclusion_field="ForeheadTimeseries",
    #                                    inclusion_value="Yes",
    #                                    personal_id_field="PersonalID",
    #                                    disturbed_field="SampleMenstruationDisturbance",
    #                                    disturbed_value="Yes")
    # print "Weighted UniFrac: %1.3f (n-disturbed samples: %d, n-undisturbed samples: %d) " % (r[1], len(r[2]), len(r[3]))
    # 
    # # <codecell>
    # 
    # r = score_ranked_adjacent_unifracs(m,
    #                                    udm,
    #                                    uh,
    #                                    inclusion_field="GutTimeseries",
    #                                    inclusion_value="Yes",
    #                                    personal_id_field="PersonalID",
    #                                    disturbed_field="SampleSicknessDisturbance",
    #                                    disturbed_value="Yes")
    # print "Unweighted UniFrac: %1.3f (n-disturbed samples: %d, n-undisturbed samples: %d) " % (r[1], len(r[2]), len(r[3]))
    # r = score_ranked_adjacent_unifracs(m,
    #                                    wdm,
    #                                    wh,
    #                                    inclusion_field="GutTimeseries",
    #                                    inclusion_value="Yes",
    #                                    personal_id_field="PersonalID",
    #                                    disturbed_field="SampleSicknessDisturbance",
    #                                    disturbed_value="Yes")
    # print "Weighted UniFrac: %1.3f (n-disturbed samples: %d, n-undisturbed samples: %d) " % (r[1], len(r[2]), len(r[3]))
    # 
    # # <codecell>
    # 
    # r = score_ranked_adjacent_unifracs(m,
    #                                    udm,
    #                                    uh,
    #                                    inclusion_field="TongueTimeseries",
    #                                    inclusion_value="Yes",
    #                                    personal_id_field="PersonalID",
    #                                    disturbed_field="SampleSicknessDisturbance",
    #                                    disturbed_value="Yes")
    # print "Unweighted UniFrac: %1.3f (n-disturbed samples: %d, n-undisturbed samples: %d) " % (r[1], len(r[2]), len(r[3]))
    # r = score_ranked_adjacent_unifracs(m,
    #                                    wdm,
    #                                    wh,
    #                                    inclusion_field="TongueTimeseries",
    #                                    inclusion_value="Yes",
    #                                    personal_id_field="PersonalID",
    #                                    disturbed_field="SampleSicknessDisturbance",
    #                                    disturbed_value="Yes")
    # print "Weighted UniFrac: %1.3f (n-disturbed samples: %d, n-undisturbed samples: %d) " % (r[1], len(r[2]), len(r[3]))
    # 
    # # <codecell>
    # 
    # r = score_ranked_adjacent_unifracs(m,
    #                                    udm,
    #                                    uh,
    #                                    inclusion_field="PalmTimeseries",
    #                                    inclusion_value="Yes",
    #                                    personal_id_field="PersonalID",
    #                                    disturbed_field="SampleSicknessDisturbance",
    #                                    disturbed_value="Yes")
    # print "Unweighted UniFrac: %1.3f (n-disturbed samples: %d, n-undisturbed samples: %d) " % (r[1], len(r[2]), len(r[3]))
    # r = score_ranked_adjacent_unifracs(m,
    #                                    wdm,
    #                                    wh,
    #                                    inclusion_field="PalmTimeseries",
    #                                    inclusion_value="Yes",
    #                                    personal_id_field="PersonalID",
    #                                    disturbed_field="SampleSicknessDisturbance",
    #                                    disturbed_value="Yes")
    # print "Weighted UniFrac: %1.3f (n-disturbed samples: %d, n-undisturbed samples: %d) " % (r[1], len(r[2]), len(r[3]))
    # 
    # # <codecell>
    # 
    # r = score_ranked_adjacent_unifracs(m,
    #                                    udm,
    #                                    uh,
    #                                    inclusion_field="ForeheadTimeseries",
    #                                    inclusion_value="Yes",
    #                                    personal_id_field="PersonalID",
    #                                    disturbed_field="SampleSicknessDisturbance",
    #                                    disturbed_value="Yes")
    # print "Unweighted UniFrac: %1.3f (n-disturbed samples: %d, n-undisturbed samples: %d) " % (r[1], len(r[2]), len(r[3]))
    # r = score_ranked_adjacent_unifracs(m,
    #                                    wdm,
    #                                    wh,
    #                                    inclusion_field="ForeheadTimeseries",
    #                                    inclusion_value="Yes",
    #                                    personal_id_field="PersonalID",
    #                                    disturbed_field="SampleSicknessDisturbance",
    #                                    disturbed_value="Yes")
    # print "Weighted UniFrac: %1.3f (n-disturbed samples: %d, n-undisturbed samples: %d) " % (r[1], len(r[2]), len(r[3]))

if __name__ == "__main__":
    main()