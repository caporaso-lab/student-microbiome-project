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

from qiime.parse import parse_mapping_file

# Globals FTW... configure at will.
schools = ['NCS', 'NAU', 'CUB']

def main():
    smp_map_f = open('smp_map.txt', 'U')
    smp_map_data, smp_map_header, smp_map_comments = \
            parse_mapping_file(smp_map_f)
    week_mapping = build_week_mapping(smp_map_data, smp_map_header)
    smp_map_f.close()

    dist_f = open('disturbance_list.txt', 'U')

    new_lines = []
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
                new_lines.append(line)
                continue

            new_line = [header]
            for cell in cells:
                if cell and '(' in cell:
                    pid, rest = cell.split('(')
                    try:
                        int(pid)
                    except:
                        raise ValueError("Encountered invalid PID '%s'." % pid)

                    if len(pid) == 1:
                        pid = '00' + pid
                    elif len(pid) == 2:
                        pid = '0' + pid

                    if len(pid) != 3:
                        raise ValueError("Encountered invalid PID '%s'." % pid)

                    # Add school to front of PID.
                    try:
                        new_pid = week_mapping[pid][0] + pid
                    except KeyError:
                        print ("Found PID '%s' in disturbance list that is "
                               "not in the SMP mapping file. Marking with two "
                               "question marks and leaving cell contents "
                               "unmodified." % pid)
                        new_line.append('??%s' % cell)
                        continue

                    weeks, rest = rest.split(';')
                    weeks = weeks.split(',')

                    valid_mapping = True
                    new_weeks = []
                    for week in weeks:
                        if week in week_mapping[pid][1]:
                            new_weeks.append(week_mapping[pid][1][week])
                        else:
                            # Mapping doesn't exist, so mark this entire cell
                            # with a question mark.
                            valid_mapping = False

                    if valid_mapping:
                        new_line.append('%s(%s;%s' %
                                        (new_pid, ','.join(new_weeks), rest))
                    else:
                        print ("Found PID '%s' in disturbance list with "
                               "week '%s' that didn't match a week in the "
                               "WeekDescription column in the SMP mapping "
                               "file. Marking with a single question mark and "
                               "leaving cell contents unmodified." %
                               (new_pid, week))
                        new_line.append('?%s(%s;%s' %
                                        (new_pid, ','.join(weeks), rest))
                else:
                    # Can't parse, so pass it through.
                    new_line.append(cell)

            new_lines.append('\t'.join(new_line))

    dist_f.close()

    new_dist_f = open('new_disturbance_list.txt', 'w')
    new_dist_f.write('\n'.join(new_lines))
    new_dist_f.close()

def build_week_mapping(map_data, map_header):
    pid_idx = map_header.index('PersonalID')
    week_desc_idx = map_header.index('WeekDescription')
    wss_idx = map_header.index('WeeksSinceStart')

    week_mapping = {}
    for row in map_data:
        school = row[pid_idx][:-3]
        pid = row[pid_idx][-3:]
        wss = row[wss_idx]

        # There are 'na' values in WeekDescription...
        try:
            week_desc = row[week_desc_idx].split('.', 1)[1]
        except:
            continue

        if school not in schools:
            continue

        if pid in week_mapping:
            if week_desc in week_mapping[pid][1]:
                if wss != week_mapping[pid][1][week_desc]:
                    raise ValueError("Encountered ambiguous mapping between "
                                     "PID '%s' and WeekDescription '%s' and "
                                     "WeeksSinceStart '%s'."
                                     % (school + pid, week_desc, wss))
            else:
                week_mapping[pid][1][week_desc] = wss
        else:
            week_mapping[pid] = (school, {week_desc: wss})

    return week_mapping


if __name__ == "__main__":
    main()
