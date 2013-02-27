#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"
__status__ = "Development"

from collections import defaultdict
from math import ceil
from qiime.format import format_mapping_file
from qiime.parse import parse_mapping_file

# Globals FTW... configure at will.
schools = ['NCS', 'NAU', 'CUB']
anti_dist_name = 'SampleAntibioticDisturbance'
sick_dist_name = 'SampleSicknessDisturbance'
menst_dist_name = 'SampleMenstruationDisturbance'

anti_headers = ['solodyn', 'clavulanate', 'difucan', 'azithromyacin',
                'bactrion', 'zipromycin', 'amoxicillin', 'doxycycline hyclate',
                'cephalexin', 'nitrofurantin', 'penicillin', 'tretinoin',
                'sulfameth', 'ciproflaxicin', 'zithromax', 'zithromyacin']

sick_headers = ['coughing', 'runny nose', 'body aches', 'fever', 'headache',
                'stress ', 'sorethroat', 'stomach ache', 'ear ache',
                'vomiting', 'diarrhea']

menst_headers = ['menstrual cycle']

def main():
    smp_map_f = open('smp_map.txt', 'U')
    smp_map_data, smp_map_header, smp_map_comments = \
            parse_mapping_file(smp_map_f)
    smp_map_f.close()

    dist_f = open('new_disturbance_list.txt', 'U')
    anti_dist_map, sick_dist_map, menst_dist_map = \
            parse_disturbance_file(dist_f)
    dist_f.close()

    pid_idx = smp_map_header.index('PersonalID')
    wss_idx = smp_map_header.index('WeeksSinceStart')

    # Add the three new columns at the end.
    new_smp_map_header = smp_map_header[:]
    new_smp_map_header.extend([anti_dist_name, sick_dist_name,
                               menst_dist_name])

    new_smp_map_data = []
    for row in smp_map_data:
        pid = row[pid_idx]
        school = pid[:-3]
        week = row[wss_idx]

        anti_dist = False
        sick_dist = False
        menst_dist = False

        # Figure out if we should try to map this sample.
        valid_sample = True

        try:
            int(pid[-3:])
        except:
            valid_sample = False

        if school not in schools:
            valid_sample = False

        try:
            week = float(week)
        except:
            valid_sample = False

        if valid_sample:
            if pid in anti_dist_map and week in anti_dist_map[pid]:
                anti_dist = True
            if pid in sick_dist_map and week in sick_dist_map[pid]:
                sick_dist = True
            if pid in menst_dist_map and week in menst_dist_map[pid]:
                menst_dist = True

        # Write out our results in three new columns.
        anti_dist_str = 'Yes' if anti_dist else 'No'
        sick_dist_str = 'Yes' if sick_dist else 'No'
        menst_dist_str = 'Yes' if menst_dist else 'No'

        new_row = row[:]
        new_row.extend([anti_dist_str, sick_dist_str, menst_dist_str])
        new_smp_map_data.append(new_row)

    new_smp_map_f = open('new_smp_map.txt', 'w')
    new_smp_map_f.write(format_mapping_file(new_smp_map_header,
                                            new_smp_map_data,
                                            smp_map_comments))
    new_smp_map_f.close()

def parse_disturbance_file(dist_f):
    anti_dist_map = defaultdict(list)
    sick_dist_map = defaultdict(list)
    menst_dist_map = defaultdict(list)

    for line in dist_f:
        line = line.strip()

        if line:
            cells = line.split('\t')
            cells = map(lambda l: l.replace('"', '').strip(), cells)
            header = cells[0]
            cells = cells[1:]

            if header == 'Disturbance type':
                # Hack to skip "detailed disturbances" worksheet header. The
                # other worksheet doesn't have a header line on its own.
                continue

            dist_map = None
            if header in anti_headers:
                dist_map = anti_dist_map
            elif header in sick_headers:
                dist_map = sick_dist_map
            elif header in menst_headers:
                dist_map = menst_dist_map

            if dist_map is not None:
                for cell in cells:
                    if cell:
                        pid, rest = cell.split('(')
                        weeks = map(int, rest.split(';')[0].split(','))
                        dist_map[pid].extend(weeks)

    return anti_dist_map, sick_dist_map, menst_dist_map


if __name__ == "__main__":
    main()
