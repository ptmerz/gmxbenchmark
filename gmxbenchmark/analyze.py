#!/usr/bin/env python3
"""
This file is part of the GROMACS benchmarking toolset.

`gmxbenchmark-analyze` reads GROMACS .log files in directory structures created
by `gmxbenchmark-prepare`, performs some simple statistical analysis of the
performance data, and presents the results in a html file.

Author: Pascal Merz <pascal.merz@me.com>

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the
  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
  Boston, MA 02110-1301 USA
"""
import argparse
import collections
import numpy as np
import os
import sys
from typing import Dict, Iterable, List, Sequence, Union


def html_head(title: str) -> str:
    """
    Returns the html head used in `gmxbenchmark-analyze`.

    Parameters
    ----------
    title: str
        The content of the `<title></title>` tag

    Returns
    -------
    str
        Head of .html file

    """
    return r"""
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<style>
table {
    font-family: arial, sans-serif;
    border-collapse: collapse;
    table-layout: auto;
}

tr {
    border-top: 1px solid #dddddd;
    border-bottom: 1px solid #dddddd;
}

td, th {
    padding: 8px;
}

th {
    background-color: #eeeeee;
}


.tooltip {
    position: relative;
    display: inline-block;
    border-bottom: 1px dotted black;
}

.tooltip .tooltiptext {
    visibility: hidden;
    background-color: black;
    color: #fff;
    text-align: center;
    border-radius: 6px;
    padding: 5px;

    /* Position the tooltip */
    position: absolute;
    z-index: 1;
    top: -5px;
    left: 105%;
}

.tooltip:hover .tooltiptext {
    visibility: visible;
}

/* Create two equal columns that floats next to each other */
.column {
  float: left;
  width: 48%;
  padding: 10px;
}

/* Clear floats after the columns */
.row:after {
  content: "";
  display: table;
  clear: both;
}
</style>
<title>""" + title + r"""</title>
</head>
<body>
"""


def html_foot() -> str:
    """
    Returns a string closing the `<body>` and `<html>` tags.

    Returns
    -------
    str
        End of .html file

    """
    return r"""

</body>
</html>
    """


def html_a(href: str, text: str, style: str = None) -> str:
    """
    Returns a html link tag.

    Parameters
    ----------
    href: str
        The url to link to
    text: str
        The link text
    style: str, optional
        CSS style string

    Returns
    -------
    str
        Link tag

    """
    link = '<a href="{:s}"'.format(href)
    if style is not None:
        link += ' style="{:s}"'.format(style)
    link += '>{:s}</a>'.format(text)
    return link


def html_header(header: str, level: int, with_id: bool = True) -> str:
    """
    Returns a html header tag.

    Parameters
    ----------
    header: str
        The header text
    level: int
        The level of the header, as in <h1> (1) or <h2> (2)
    with_id: bool, optional
        Whether a `id=` tag is added. If yes, use `header` as id,
        replacing spaces by underscores

    Returns
    -------
    str
        Header tag

    """
    assert level >= 1
    if with_id:
        anchor = header.replace(' ', '_')
        return '<h{0:d} id="{2:s}">{1:s}</h{0:d}>\n'.format(level, header, anchor)
    else:
        return '<h{0:d}>{1:s}</h{0:d}>\n'.format(level, header)


def html_paragraph(string: str) -> str:
    """
    Return the argument string surrounded by `<p></p>` tags.

    Parameters
    ----------
    string: str
        The content of the paragraph

    Returns
    -------
    str
        Paragraph tag

    """
    return '<p>{:s}</p>\n'.format(string)


def html_codeblock(code: str) -> str:
    """
    Return a html block displaying a pre-formatted code block.

    Parameters
    ----------
    code: str
        Formatted code block

    Returns
    -------
    str
        Code tag

    """
    return '<pre><code>\n' \
           '{:s}\n' \
           '</code></pre>'.format(code)


def html_table(
        header: List[str], rows: List[List[str]], footer: List[str] = None,
        alternate: bool = False, header_spans2: bool = False) -> str:
    """
    Returns a html table.

    Parameters
    ----------
    header: List[str]
        The header row
    rows: List[List[str]]
        A list of table body rows
    footer: List[str], optional
        The footer row, similar in appearance to the header row
    alternate: bool, optional
        Whether the columns should alternate in left and right alignment,
        useful for alternating value - error estimate pairs
    header_spans2: bool, optional
        Whether the header row (except the first column) should span two
        columns per cell.

    Returns
    -------
    str
        The table

    """
    table = ['<table>']
    if header:
        table.append(' ' * 2 + '<tr style="text-align: left;">')
        for idx, h in enumerate(header):
            add = ''
            if header_spans2 and idx > 0:
                add = ' colspan="2"'
            table.append(' ' * 4 + '<th{:s}>'.format(add) + h + '</th>')
        table.append(' ' * 2 + '</tr>')
    for row in rows:
        table.append(' ' * 2 + '<tr>')
        if alternate:
            align = ' style="text-align: left; padding-left: 2px;"'
            other = ' style="text-align: right; padding-right: 2px;"'
        else:
            align = ''
            other = ''
        for r in row:
            table.append(' ' * 4 + '<td{:s}>'.format(align) + str(r) + '</td>')
            align, other = other, align
        table.append(' ' * 2 + '</tr>')
    if footer:
        table.append(' ' * 2 + '<tr>')
        if alternate:
            align = ' text-align: left; padding-left: 2px;'
            other = ' text-align: right; padding-right: 2px;'
        else:
            align = ''
            other = ''
        for f in footer:
            table.append(' ' * 4 + '<td style="background-color: #eeeeee; font-weight: bold;{:s}">'.format(align) +
                         f + '</td>')
            align, other = other, align
        table.append(' ' * 2 + '</tr>')
    table.append('</table>')

    return '\n'.join(table)


def html_two_cols(entry1: str, entry2: str) -> str:
    """
    Returns a html construct allowing to have two columns.

    Parameters
    ----------
    entry1: str
        Content of column 1
    entry2: str
        Content of column 2

    Returns
    -------
    str
        The double column html construct

    """
    assert isinstance(entry1, str) and isinstance(entry2, str)
    html_list = list()
    html_list.append('<div class="row">')
    html_list.append('  <div class="column">')
    html_list.append(' '*4 + entry1)
    html_list.append('  </div>')
    html_list.append('  <div class="column">')
    html_list.append(' '*4 + entry2)
    html_list.append('  </div>')
    html_list.append('</div>')
    return '\n'.join(html_list)


class GroFile:
    """Objects of this class contain information from .gro files

    Objects of this class read a GROMACS .gro file, and make available the
    relevant information - currently only the number of atoms and the box
    vector.
    """
    def __init__(self, grofile: Iterable[str]) -> None:
        """
        Initialize the GroFile object.

        Parameters
        ----------
        grofile: Iterable[str]
            File handler or content of GROMACS .gro file

        """
        lines = [line.rstrip('\n') for line in grofile if line.strip()]
        self._natoms = int(lines[1])
        self._box = np.array([float(n) for n in lines[-1].split()[:3]])

    def get_natoms(self) -> int:
        """
        Return the number of atoms found in the .gro file.

        Returns
        -------
        int
            Number of atoms

        """
        return self._natoms

    def get_box(self) -> np.ndarray:
        """
        Return the box vector found in the .gro file

        Returns
        -------
        np.ndarray(shape=(1, 3), dtype=float)
            The box vector

        """
        return self._box


class LogFile:
    """Information from .log files

    Objects of this class read a GROMACS .log file, and make available the
    relevant information - currently the different performance measurements
    and the version and the path of the executable.
    """
    def __init__(self, logfile: Iterable[str]) -> None:
        """
        Initialize the LogFile object. This reads the file and extracts the
        relevant information.

        Parameters
        ----------
        logfile: Iterable[str]
            File handler or content of GROMACS .log file

        """
        self._performanceblock = list()
        self._exe = ''
        self._version = ''
        footer_found = False
        for line in logfile:
            if line.startswith('Executable:'):
                self._exe = line.split()[-1]
            if line.startswith('GROMACS version:'):
                self._version = line.split()[-1]
            if 'M E G A - F L O P S   A C C O U N T I N G' not in line and not footer_found:
                continue
            footer_found = True
            self._performanceblock.append(line)

        self._flops_block = dict()
        self._read_flops_block()
        
        self._time_block = dict()
        self._read_time_block()

        self._walltime = 0
        self._coretime = 0
        self._read_times()

    def get_flops_block(self) -> Dict[str,
                                      Union[List[str],
                                            Dict[str, Union[List[str],
                                                            Dict[str,
                                                                 np.ndarray],
                                                            np.ndarray]]]]:
        """
        Return block representing the FLOPS accounting data.

        Returns
        -------
            Representation of FLOPS accounting data
        """
        return self._flops_block

    def get_time_block(self) -> Dict[str,
                                     Union[List[str],
                                           Dict[str, Union[List[str],
                                                           Dict[str, np.ndarray],
                                                           np.ndarray]]]]:
        """
        Return block representing the time accounting data.

        Returns
        -------
            Representation of time accounting data
        """
        return self._time_block

    def get_performance_block(self) -> List[str]:
        """
        Return the entire performance log block (tail of GROMACS .log file)

        Returns
        -------
        List[str]
            The performance block

        """
        return self._performanceblock

    def get_version(self) -> str:
        """
        Return the version string found in the .log file.

        Returns
        -------
        str
            The GROMACS version

        """
        return self._version

    def get_exe(self) -> str:
        """
        Return the path to the GROMACS executable from the .log file.

        Returns
        -------
        str
            The path to the GROMACS executable

        """
        return self._exe

    def get_walltime(self) -> float:
        """
        Return the wall time found in the log file.

        Returns
        -------
        float
            The wall time found in the log file

        """
        return self._walltime

    def get_coretime(self) -> float:
        """
        Return the core time found in the log file.

        Returns
        -------
        float
            The core time found in the log file

        """
        return self._coretime

    def _read_flops_block(self) -> None:
        """
        Extracts and stores the FLOPS accounting from the performance block.
        """
        separator = '-----------------------------------------------------------------------------'
        self._flops_block['title'] = 'MEGA-FLOPS ACCOUNTING'
        seps = [idx for idx, line in enumerate(self._performanceblock) if separator in line]
        self._flops_block['desc'] = [line for line in self._performanceblock[1:seps[0] - 1] if line.strip()]
        self._flops_block['table'] = dict()

        self._flops_block['table']['header'] = ['Computing:', 'M-Number', 'M-Flops', '% Flops']
        # A little sanity check
        if ' {:32s} {:>16s} {:>15s}  {:>7s}\n'.format(
                *self._flops_block['table']['header']) != self._performanceblock[seps[0] - 1]:
            raise RuntimeError('Log file format not recognized.')

        self._flops_block['table']['body'] = collections.OrderedDict()
        for line in self._performanceblock[seps[0] + 1:seps[1]]:
            line = line.split()
            key = ' '.join(line[:-3])
            self._flops_block['table']['body'][key] = np.array(
                [float(line[-3]), float(line[-2]), float(line[-1])])
        line = self._performanceblock[seps[1] + 1].split()
        self._flops_block['table']['footer'] = np.array([float(line[-2]), float(line[-1])])

    def _read_time_block(self) -> None:
        """
        Extracts and stores the time accounting from the performance block.
        """
        separator = '-----------------------------------------------------------------------------'
        self._time_block['title'] = 'REAL CYCLE AND TIME ACCOUNTING'
        title_idx = next(idx for idx, line in enumerate(self._performanceblock)
                         if 'R E A L   C Y C L E   A N D   T I M E   A C C O U N T I N G' in line)
        seps = [idx + title_idx for idx, line in enumerate(self._performanceblock[title_idx:]) if separator in line]
        self._time_block['desc'] = [line for line in self._performanceblock[title_idx + 1:seps[0] - 2] if line.strip()]
        self._time_block['table'] = dict()

        self._time_block['table']['header'] = ['Computing:', 'Num Ranks', 'Num Threads', 'Call Count',
                                               'Wall time (s)', 'Giga-Cycles total sum', 'Giga-Cycles %']
        # A little sanity check
        if ((' Computing:          Num   Num      Call    Wall time         Giga-Cycles\n'
             != self._performanceblock[seps[0] - 2]) or
            ('                     Ranks Threads  Count      (s)         total sum    %\n'
             != self._performanceblock[seps[0] - 1])):
            raise RuntimeError('Log file format not recognized.')

        self._time_block['table']['body'] = collections.OrderedDict()
        for line in self._performanceblock[seps[0] + 1:seps[1]]:
            key = line[:20].strip()
            line = line.split()
            try:
                self._time_block['table']['body'][key] = np.array(
                    [int(line[-6]), int(line[-5]), int(line[-4]),
                     float(line[-3]), float(line[-2]), float(line[-1])])
            except IndexError or ValueError:
                self._time_block['table']['body'][key] = np.array(
                    [-1, -1, -1,
                     float(line[-3]), float(line[-2]), float(line[-1])])
        line = self._performanceblock[seps[1] + 1].split()
        self._time_block['table']['footer'] = np.array([float(line[-3]), float(line[-2]), float(line[-1])])

    def _read_times(self) -> None:
        """
        Extracts and stores the wall time and the core time from the performance block.
        """
        for line in self._performanceblock:
            if line.startswith('       Time:'):
                self._coretime = float(line.split()[1])
                self._walltime = float(line.split()[2])


class TimingStatistics:
    """Statistics over informations from .log files

    Objects of this class evaluate the performance informations from multiple `LogFile`
    objects to perform simple statistical analysis of the performance of a simulation
    over multiple repetitions. It also checks that all log-files were written by the
    same GROMACS executable, and offers access to the GROMACS version and executable
    path used.
    """
    def __init__(self, logfiles: Sequence[LogFile]) -> None:
        """
        Initialize TimingStatistics object from a list of LogFile objects. This function
        checks that all version strings and executable paths are identical and stores this
        information. Then, it proceeds to save the information of the flops and time
        accounting blocks of the LogFile objects as numpy arrays, ready for further
        statistical analysis.

        Parameters
        ----------
        logfiles: Sequence[LogFile]
            List of LogFile objects

        """
        versions = [log.get_version() for log in logfiles]
        executables = [log.get_exe() for log in logfiles]

        if not len(set(versions)) == 1:
            raise RuntimeError('Log files were not produced with the same GROMACS version')
        if not len(set(executables)) == 1:
            raise RuntimeError('Log files were not produced with the same GROMACS executable')

        self._version = versions[0]
        self._executable = executables[0]

        self._performanceblocks = [log.get_performance_block() for log in logfiles]
        
        self._walltimes = np.array([log.get_walltime() for log in logfiles])
        self._coretimes = np.array([log.get_coretime() for log in logfiles])

        # Collect flops
        self._flops = dict()
        for field in ['title', 'desc']:
            self._flops[field] = logfiles[0].get_flops_block()[field]
        self._flops['table'] = dict()
        self._flops['table']['header'] = logfiles[0].get_flops_block()['table']['header']
        self._flops['table']['body'] = collections.OrderedDict()
        
        for key in logfiles[0].get_flops_block()['table']['body']:
            self._flops['table']['body'][key] = np.zeros((len(logfiles),
                                                          len(logfiles[0].get_flops_block()['table']['body'][key])))
        for idx, timing in enumerate(logfiles):
            for key in timing.get_flops_block()['table']['body']:
                assert key in self._flops['table']['body']
                self._flops['table']['body'][key][idx, :] = timing.get_flops_block()['table']['body'][key]

        self._flops['table']['footer'] = np.array([timing.get_flops_block()['table']['footer'] for timing in logfiles])

        # Collect time
        self._time = dict()
        for field in ['title', 'desc']:
            self._time[field] = logfiles[0].get_time_block()[field]
        self._time['table'] = dict()
        self._time['table']['header'] = logfiles[0].get_time_block()['table']['header']
        self._time['table']['body'] = collections.OrderedDict()
        for key in logfiles[0].get_time_block()['table']['body']:
            self._time['table']['body'][key] = np.zeros((len(logfiles),
                                                         len(logfiles[0].get_time_block()['table']['body'][key])))
        for idx, timing in enumerate(logfiles):
            for key in timing.get_time_block()['table']['body']:
                assert key in self._time['table']['body']
                self._time['table']['body'][key][idx, :] = timing.get_time_block()['table']['body'][key]

        self._time['table']['footer'] = np.array([timing.get_time_block()['table']['footer'] for timing in logfiles])

    def get_flopsblock(self) -> Dict[str,
                                     Union[str,
                                           List[str],
                                           Dict[str,
                                                Union[List[str], List[List[str]]]]]]:
        """
        Return the FLOPS block ready for output in a html table. Includes the title and description
        of the block, plus a table including header, body and footer. The table values have been
        calculated as arithmetic mean of the values extracted from the single log files, and are
        reported with the standard deviation as error estimate.

        Returns
        -------
            Dictionary containing title, description, and result table

        """
        flops = dict()
        flops['title'] = self._flops['title']
        flops['desc'] = self._flops['desc']
        
        header = self._flops['table']['header']
        rows = []
        for key in self._flops['table']['body']:
            line = np.array([
                [mean, var] for mean, var in zip(
                    self._flops['table']['body'][key].mean(axis=0),
                    self._flops['table']['body'][key].std(axis=0))]).flatten()
            rows.append([
                key,
                '{:.6f}'.format(line[0]), '&plusmn;{:.6f}'.format(line[1]).rstrip('0').rstrip('.'),
                '{:.3f}'.format(line[2]), '&plusmn;{:.3f}'.format(line[3]).rstrip('0').rstrip('.'),
                '{:.1f}'.format(line[4]), '&plusmn;{:.1f}'.format(line[5]).rstrip('0').rstrip('.')])

        footer = np.array([
                [mean, var] for mean, var in zip(
                    self._flops['table']['footer'].mean(axis=0),
                    self._flops['table']['footer'].std(axis=0))]).flatten()
        footer = [
            'Total', '', '',
            '{:.3f}'.format(footer[0]), '&plusmn;{:.3f}'.format(footer[1]).rstrip('0').rstrip('.'),
            '{:.1f}'.format(footer[2]), '&plusmn;{:.1f}'.format(footer[3]).rstrip('0').rstrip('.')]
        
        flops['table'] = html_table(header, rows, footer, True, True)
        return flops

    def get_timeblock(self) -> Dict[str,
                                    Union[str,
                                          List[str],
                                          Dict[str,
                                               Union[List[str], List[List[str]]]]]]:
        """
        Return the time block ready for output in a html table. Includes the title and description
        of the block, plus a table including header, body and footer. The table values have been
        calculated as arithmetic mean of the values extracted from the single log files, and are
        reported with the standard deviation as error estimate.

        Returns
        -------
            Dictionary containing title, description, and result table

        """
        time = dict()
        time['title'] = self._time['title']
        time['desc'] = self._time['desc']
        
        header = self._time['table']['header']
        rows = []
        for key in self._time['table']['body']:
            line = np.array([
                [mean, var] for mean, var in zip(
                    self._time['table']['body'][key].mean(axis=0),
                    self._time['table']['body'][key].std(axis=0))]).flatten()
            row = [
                key,
                '{:.1g}'.format(line[0]), '&plusmn;{:.1f}'.format(line[1]).rstrip('0').rstrip('.'),
                '{:.1g}'.format(line[2]), '&plusmn;{:.1f}'.format(line[3]).rstrip('0').rstrip('.'),
                '{:.0f}'.format(line[4]), '&plusmn;{:.0f}'.format(line[5]),
                '{:.3f}'.format(line[6]), '&plusmn;{:.3f}'.format(line[7]),
                '{:.3f}'.format(line[8]), '&plusmn;{:.3f}'.format(line[9]),
                '{:.1f}'.format(line[10]), '&plusmn;{:.1f}'.format(line[11])
            ]
            if row[1] == '-1':
                row[1] = ''
                row[2] = ''
            if row[3] == '-1':
                row[3] = ''
                row[4] = ''
            if row[5] == '-1':
                row[5] = ''
                row[6] = ''
            rows.append(row)
        
        footer = np.array([
                [mean, var] for mean, var in zip(
                    self._time['table']['footer'].mean(axis=0),
                    self._time['table']['footer'].std(axis=0))]).flatten()
        footer = [
            'Total', '', '', '', '', '', '',
            '{:.3f}'.format(footer[0]), '&plusmn;{:.3f}'.format(footer[1]),
            '{:.3f}'.format(footer[2]), '&plusmn;{:.3f}'.format(footer[3]),
            '{:.1f}'.format(footer[4]), '&plusmn;{:.1f}'.format(footer[5])
        ]

        time['table'] = html_table(header, rows, footer, True, True)
        return time

    def get_wall_time(self) -> np.ndarray:
        """
        Return the arithmetic mean of the wall time read in the single log files,
        as well as the standard deviation as error estimate.

        Returns
        -------
        np.ndarray(shape=(1,2), dtype=float)
            Arithmetic mean and standard deviation of the wall time.

        """
        return np.array((self._walltimes.mean(), self._walltimes.std()))

    def get_core_time(self) -> np.ndarray:
        """
        Return the arithmetic mean of the core time read in the single log files,
        as well as the standard deviation as error estimate.

        Returns
        -------
        np.ndarray(shape=(1,2), dtype=float)
            Arithmetic mean and standard deviation of the core time.

        """
        return np.array((self._coretimes.mean(), self._coretimes.std()))

    def get_performance_blocks(self) -> List[str]:
        """
        Return a list of formatted strings containing the performance blocks of
        all .log files.

        Returns
        -------
        List[str]
            List of performance blocks

        """
        return [''.join(t) for t in self._performanceblocks]

    def get_version(self) -> str:
        """
        Return the version string, which was checked to be identical for
        all .log files.

        Returns
        -------
        str
            Version

        """
        return self._version

    def get_exe(self) -> str:
        """
        Return the executable path, which was checked to be identical for
        all .log files.

        Returns
        -------
        str
            Executable path

        """
        return self._executable


def write_dropdown_diff(diff_options: List[Union[str, List[str]]],
                        common_options: List[Union[str, List[str]]]) -> str:
    """
    Create a html code block which allows to specify two anchors via dropdown
    lists, and replace the content of two <div> tags with this the content
    of <div> tags matching the chosen anchors. This implements an easy way
    to compare timing information between different runs.

    Parameters
    ----------
    diff_options: List[Union[str, List[str]]]
        Options for one or more dropdown lists which are different for the
        two targets
    common_options: List[Union[str, List[str]]]
        Options for one or more dropdown lists which are identical for the
        two targets

    Returns
    -------
    str
        The html code block

    """
    if not isinstance(diff_options[0], list):
        diff_options = [diff_options]
    if not isinstance(common_options[0], list):
        common_options = [common_options]

    dropdowns_diff = [
        ['<option value="{0:s}">{0:s}</option>'.format(option)
         for option in diff_option]
        for diff_option in diff_options
    ]
    dropdowns_common = [
        ['<option value="{0:s}">{0:s}</option>'.format(option)
         for option in common_option]
        for common_option in common_options
    ]

    common = list()
    for ndd, dropdown in enumerate(dropdowns_common):
        common.append('<select id="dd-common-{:d}">'.format(ndd))
        common.append('\n'.join(dropdown))
        common.append('</select>')
    common.append('<button onclick="showDiff()">Compare</button>')
    common.append('<button onclick="clearDiff()">Clear</button>')

    diff1 = list()
    for ndd, dropdown in enumerate(dropdowns_diff):
        diff1.append('<select id="dd-diff1-{:d}">'.format(ndd))
        diff1.append('\n'.join(dropdown))
        diff1.append('</select>')
    diff2 = list()
    for ndd, dropdown in enumerate(dropdowns_diff):
        diff2.append('<select id="dd-diff2-{:d}">'.format(ndd))
        diff2.append('\n'.join(dropdown))
        diff2.append('</select>')

    script = list()
    script.append('<script>')
    script.append('function showDiff() {')

    script.append('  var div_diff1 = document.getElementById("diff1");')
    script.append('  var div_diff2 = document.getElementById("diff2");')
    script.append('  var common_id = "";')
    script.append('  var diff1_id = "";')
    script.append('  var diff2_id = "";')

    script.append('  for (var idx = 0; idx < {:d}; idx++) {{'.format(len(dropdowns_common)))
    script.append('    var dropdown = document.getElementById("dd-common-" + String(idx));')
    script.append('    common_id += "-" + dropdown.value;')
    script.append('  }')
    script.append('  for (var idx = 0; idx < {:d}; idx++) {{'.format(len(dropdowns_diff)))
    script.append('    var dropdown = document.getElementById("dd-diff1-" + String(idx));')
    script.append('    diff1_id += "-" + dropdown.value;')
    script.append('  }')
    script.append('  for (var idx = 0; idx < {:d}; idx++) {{'.format(len(dropdowns_diff)))
    script.append('    var dropdown = document.getElementById("dd-diff2-" + String(idx));')
    script.append('    diff2_id += "-" + dropdown.value;')
    script.append('  }')

    script.append('  var source_diff1 = document.getElementById("timing" + diff1_id + common_id);')
    script.append('  var source_diff2 = document.getElementById("timing" + diff2_id + common_id);')
    script.append('  div_diff1.innerHTML = source_diff1.innerHTML;')
    script.append('  div_diff2.innerHTML = source_diff2.innerHTML;')
    script.append('}')

    script.append('function clearDiff() {')
    script.append('  var div_diff1 = document.getElementById("diff1");')
    script.append('  var div_diff2 = document.getElementById("diff2");')
    script.append('  div_diff1.innerHTML = "Press \\"Compare\\" to start diff."')
    script.append('  div_diff2.innerHTML = "Press \\"Compare\\" to start diff."')
    script.append('}')
    script.append('</script>')

    diff_html = list()
    diff_html.append(
        html_paragraph('\n'.join(common))
    )

    diff_html.append(html_two_cols(
        html_paragraph('\n'.join(diff1)) +
        '\n<div id="diff1">\nPress "Compare" to start diff.\n</div>',
        html_paragraph('\n'.join(diff2)) +
        '\n<div id="diff2">\nPress "Compare" to start diff.\n</div>'
    ))

    diff_html.append('\n'.join(script))

    return '\n'.join(diff_html)


def analyze(args: Sequence[str]) -> int:
    """
    `analyze` reads GROMACS .log files in directory structures created
    by `gmxbenchmark-prepare`, performs some simple statistical analysis of the
    performance data, and presents the results in a html file.

    Parameters
    ----------
    args: Sequence[str]
        Command line arguments to be parsed by argparse

    Returns
    -------
    int

    """
    parser = argparse.ArgumentParser(
        description='Analyze benchmarks for GROMACS',
        prog='gmxbenchmark-prepare',
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('directories', type=str, metavar='dir', nargs='+',
                        help='One or more run directories.')
    parser.add_argument('--nsystemsizes', type=int, metavar='N', default=1,
                        help='How many system sizes were simulated. Default: 1.')
    parser.add_argument('--run_args', type=str, metavar='"\\-gmx_arg X"', action='append',
                        help='Every use of --run_args denotes an independent run of `gmx mdrun` '
                             'with the given\narguments. Note: Prepend a \'\\\'.')
    parser.add_argument('--nreps', type=int, metavar='N', default=5,
                        help='The number of times the benchmarking runs are repeated, with the timings being\n'
                             'averaged over all repetitions. Default: 5.')
    parser.add_argument('--details', action='store_true',
                        help='Add more details to the resulting html file.')

    args = parser.parse_args(args)

    run_args = ['']
    if args.run_args is not None:
        run_args = args.run_args

    to_plot = []        # [ [np.vector, desc (dX, aY)] ]
    tables = []         # [ [title, [title-str, title-str, ...], [str, str, ...] ] ]
    detail_tables = []  # [ [ anchor, dict, dict ] ]
    code_blocks = []    # [ [ anchor, [ str ] ] ]
    natoms = dict()
    boxes = dict()
    executable = dict()
    version = dict()

    for ndir, directory in enumerate(args.directories):
        tables.append(['Executable: ' + directory])
        tableheader = [''] + ['Size ' + str(i+1) for i in range(args.nsystemsizes)]
        tables[-1].append(tableheader)
        for nra, ra in enumerate(run_args):
            code = 'exe{:d}-arg{:d}'.format(ndir + 1, nra + 1)
            timing_row = ['<b>{:s}</b>'.format(code)]
            plot_row = list()
            error_row = list()
            for nsystemsize in range(args.nsystemsizes):
                logfiles = []
                for nrep in range(args.nreps):
                    logfiles.append(
                        LogFile(open(os.path.join(directory, 'benchmark_{:d}/run_{:d}_{:d}.log'.format(
                            nsystemsize + 1, nra + 1, nrep + 1)))))
                timing_stats = TimingStatistics(logfiles)

                grofile = GroFile(
                    open(os.path.join(directory, 'benchmark_{:d}/run_{:d}_1.gro'.format(nsystemsize + 1, nra + 1))))

                normalized_time = timing_stats.get_core_time() / grofile.get_natoms() * 1000

                anchor = code + '-sys{:d}'.format(nsystemsize + 1)
                detail_tables.append([anchor, timing_stats.get_flopsblock(), timing_stats.get_timeblock()])
                code_blocks.append([anchor, timing_stats.get_performance_blocks()])
                timing_row.append('{:.2f}'.format(normalized_time[0]))
                timing_row.append('&plusmn;{:.2f}&nbsp;'.format(normalized_time[1]) +
                                  html_a('#' + anchor, '&#9432;', style='text-decoration: none;'))
                plot_row.append(normalized_time[0])
                error_row.append(normalized_time[1])
                if nsystemsize not in natoms:
                    natoms[nsystemsize] = grofile.get_natoms()
                elif natoms[nsystemsize] != grofile.get_natoms():
                    raise RuntimeError('Results of same system size differ in number of atoms')
                if nsystemsize not in boxes:
                    boxes[nsystemsize] = np.array([grofile.get_box()])
                else:
                    boxes[nsystemsize] = np.append(boxes[nsystemsize], [grofile.get_box()], axis=0)
                if directory not in executable:
                    executable[directory] = timing_stats.get_exe()
                elif executable[directory] != timing_stats.get_exe():
                    raise RuntimeError('Results in same folder were not produced by same GROMACS executable')
                if directory not in version:
                    version[directory] = timing_stats.get_version()
                elif version[directory] != timing_stats.get_version():
                    raise RuntimeError('Results in same folder were not produced by same GROMACS version')

            to_plot.append((np.array(plot_row), np.array(error_row), code))
            tables[-1].append(timing_row)

    do_plot = args.nsystemsizes > 2
    plot_name = None
    if do_plot:
        import matplotlib.pyplot as plt
        plot_name = 'benchmarks.png'
        fig, ax = plt.subplots()
        for line in to_plot:
            ax.errorbar(range(1, len(line[0])+1), line[0], yerr=line[1], label=line[2])
        ax.set(xlabel='System size', ylabel='Core time / natoms [ms]')
        ax.legend()
        fig.savefig(plot_name, dpi=300)

    with open('analysis.html', 'w') as html_file:
        html_file.write(html_head('GROMACS benchmark'))
        html_file.write(html_header('GROMACS benchmark', 1))
        html_file.write(html_header('Key', 2))
        html_file.write(html_paragraph('Each combination (executable - arguments - system size) was ran '
                                       '{:d} times to improve statistics.'.format(args.nreps)))
        col1 = list()
        col1.append(html_header('Executable', 3))
        col1.append(html_table(['Key', 'Name', 'Version', 'Full path'],
                               [['<b>exe{:d}</b>'.format(ndir + 1), directory,
                                 version[directory], executable[directory]]
                                for ndir, directory in enumerate(args.directories)]))
        col1.append(html_header('Arguments', 3))
        col1.append(html_table(['Key', 'Run arguments'],
                               [['<b>arg{:d}</b>'.format(nra + 1), ra[1:]] for nra, ra in enumerate(run_args)]))
        col2 = list()
        col2.append(html_header('System sizes', 3))
        col2.append(html_table(['Key', 'Size of system'],
                               [['<b>sys{:d}</b>'.format(idx+1),
                                 '# atoms: {:d}, average box size = [{:g}, {:g}, {:g}]'.format(
                                     natoms[idx], *boxes[idx].mean(axis=0))]
                                for idx in range(args.nsystemsizes)]))
        col2.append(html_paragraph('Note that the average box size is taken over the final configuration of the runs, '
                                   'not over the entire runs. For simulations in which the box is not held constant '
                                   'this is not an accurate measure of the box size over the '
                                   'simulation, but only a help to estimate the approximate system size used for '
                                   'benchmarking.'))
        html_file.write(html_two_cols('\n'.join(col1), '\n'.join(col2)))

        html_file.write(html_header('Plot', 2))
        if not do_plot:
            html_file.write(html_paragraph('No plot produced (less than 3 system sizes)'))
        else:
            html_file.write(html_a(plot_name,
                                   '<img src="{:s}" alt="Core time plot" width="600" />\n'.format(plot_name)))

        html_file.write(html_header('Core times', 2))
        for table in tables:
            html_file.write(html_header(table[0], 3))
            html_file.write(html_header('Core time per atom [ms]', 4, False))
            html_file.write(html_table(table[1], table[2:], alternate=True, header_spans2=True))

        html_file.write(html_header('Compare timings', 2))
        html_file.write(write_dropdown_diff(
            ['exe{:d}'.format(ndir+1) for ndir in range(len(args.directories))],
            [['arg{:d}'.format(narg+1) for narg in range(len(run_args))],
             ['sys{:d}'.format(nsys+1) for nsys in range(args.nsystemsizes)]]))

        html_file.write(html_header('Timings', 2))
        for table in detail_tables:
            anchor = table[0]
            html_file.write('<div id=timing-{:s}>'.format(anchor))
            html_file.write(html_header(anchor, 3))
            html_file.write(html_a('#Core_times', 'Back to tables'))
            for idx in range(1, 3):
                html_file.write(html_header(table[idx]['title'], 4, False))
                html_file.write(html_codeblock('\n'.join(table[idx]['desc'])))
                html_file.write(table[idx]['table'])
            if args.details:
                html_file.write(html_paragraph('<b>Raw output:</b> ' + ' | '.join(
                    [html_a('#{:s}-run{:d}'.format(anchor, run+1), '{:s}-run{:d}'.format(anchor, run+1))
                     for run in range(args.nreps)]
                )))
            html_file.write('</div>')

        if args.details:
            html_file.write(html_header('Raw output', 2))
            for cb in code_blocks:
                anchor = cb[0]
                for run, codeblock in enumerate(cb[1]):
                    html_file.write(html_header(cb[0] + '-run{:d}'.format(run+1), 3))
                    html_file.write(html_a('#' + anchor, 'Back to timings'))
                    html_file.write(html_codeblock(codeblock))

        html_file.write(html_foot())

    return 0


def main():
    sys.exit(analyze(sys.argv[1:]))


if __name__ == '__main__':
    main()
