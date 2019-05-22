#!/usr/bin/env python3
"""
This file is part of the GROMACS benchmarking toolset.

`gmxbenchmark-prepare` prepares files and scripts for equilibration, run and
analysis of simulations to determine the relative performance of different
GROMACS excutables.

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
import pathlib  # require py3.5
import shutil
import sys
from typing import Dict, List


def read_mdp(mdp_file_name: str, mdp_dict: Dict[str, str] = None) -> Dict[str, str]:
    """
    Read a GROMACS .mdp file, return a (ordered) dict containing the
    keys and values read in the mdp file.

    Parameters
    ----------
    mdp_file_name: str
        File name of the .mdp file
    mdp_dict: dict, optional
        Fill mdp entries in this dict, overwrite values if already existing

    Returns
    -------
    mdp_dict: dict
        Dictionary containing keys and values read in the mdp file

    """
    if mdp_dict is None:
        mdp_dict = collections.OrderedDict()
    with open(mdp_file_name) as mdp_file:
        for line in mdp_file:
            key, value = line.split('=', 1)
            mdp_dict[key.strip()] = value.strip()
    return mdp_dict


def write_mdp(mdp_file_name: str, mdp_dict: Dict[str, str]) -> None:
    """
    Writes entries of a dictionary (key, value) to a GROMACS .mdp file

    Parameters
    ----------
    mdp_file_name: str
        File name of the .mdp file
    mdp_dict: Dict[str, str]
        mdp entries to be written

    """
    with open(mdp_file_name, 'w') as mdp_file:
        for key in mdp_dict:
            mdp_file.write('{:24s} = {:s}'.format(key, mdp_dict[key]))


def get_box_size(gro_file_name: str) -> np.ndarray(dtype=float, shape=(1,3)):
    """
    Returns the box size in a GROMACS .gro file

    Parameters
    ----------
    gro_file_name: str
        File name of the .gro file

    Returns
    -------
    box: numpy.ndarray(dtype=float, shape=(1,3))
        The box size found in the .gro file

    """
    with open(gro_file_name) as gro_file:
        lines = [line.strip() for line in gro_file if line.strip()]
    return np.array([float(n) for n in lines[-1].split()[:3]])


def strip_topology(read_file: str, write_file: str) -> None:
    """
    Reads GROMACS .top file from one location, and writes it to another location
    leaving the `[ molecules ]` block empty. This is useful in connection
    with the `gmx solvate` command, which creates a box of prescribed
    dimensions and appends a line containing the number of molecules in
    the created system to the topology file.

    This function prints a warning if more than one molecule type is
    detected in the topology, as the combination of several molecules,
    `gmx solvate` and the benchmarking tool has not been tested.

    Parameters
    ----------
    read_file: str
        File name of the .top file to read from
    write_file
        File name of the .top file to write to

    """
    lines = []
    molecules = []
    with open(read_file) as read_top:
        molecules_block = False
        for line in read_top:
            if '[ molecules ]' not in line and not molecules_block:
                lines.append(line)
                continue
            if '[ molecules ]' in line:
                molecules_block = True
                lines.append(line)
                continue
            if line.startswith(';') or not line.strip():
                lines.append(line)
                continue

            molecules.append(line)

    if len(molecules) > 1:
        print('WARNING: More than one type of molecule detected - this has not been tested.\n')

    with open(write_file, 'w') as write_top:
        for line in lines:
            write_top.write(line)


class Script:
    """Helper class to write the different scripts

    Objects of this class offer methods to set various options about the
    system to be simulated and the command line options for the equilibration
    and simulation runs. They also offer a `write_scripts` function which
    write equilibration, run and analysis scripts to a user-defined
    directory.
    """
    def __init__(self) -> None:
        """
        Construct empty object
        """
        self._directory = None
        self._gmx_exe = None
        self._gmx_names = None
        self._nreps = None
        self._multiply = None
        self._num_system_sizes = None
        self._min_box = None
        self._system_name = None
        self._systems_set = False
        self._mdrun_args = ''
        self._eq_args = ''
        self._run_args = ['']
        self._mdp_path = ''

    def set_directory(self, directory: str) -> None:
        """
        Set the directory the scripts will get written to. Calling
        this function before using the `write_scripts` function is
        mandatory.

        Parameters
        ----------
        directory: str
            The directory, needs to be existing

        """
        assert os.path.exists(directory)
        self._directory = os.path.abspath(directory)

    def set_gmx(self, gmx_exe: List[str], gmx_names: List[str]) -> None:
        """
        Set the GROMACS executables used in the scripts, and a list
        containing there names (used for directory naming and output).
        The two lists are required to have the same length. Calling this
        function before using the `write_scripts` function is mandatory.

        Parameters
        ----------
        gmx_exe: List[str]
            List of GROMACS executables used in the scripts
        gmx_names: List[str]
            List of names for the executables in `gmx_exe`

        """
        assert len(gmx_exe) == len(gmx_names)
        for idx in range(len(gmx_exe)):
            if os.path.exists(gmx_exe[idx]):
                gmx_exe[idx] = os.path.abspath(gmx_exe[idx])
        self._gmx_exe = gmx_exe
        self._gmx_names = gmx_names

    def set_mdrun_args(self, args: str) -> None:
        """
        Set command line arguments used for every call to `gmx mdrun`
        during benchmarking runs. Calling this function before using
        `write_scripts` is optional.

        Parameters
        ----------
        args: str
            A string containing any number of command line arguments to `gmx mdrun`

        """
        self._mdrun_args = args

    def set_eq_args(self, args: str) -> None:
        """
        Set command line arguments used for every call to `gmx mdrun`
        during equilibration runs. Calling this function before using
        `write_scripts` is optional.

        Parameters
        ----------
        args: str
            A string containing any number of command line arguments to `gmx mdrun`

        """
        self._eq_args = args

    def set_run_args(self, args_list: List[str]) -> None:
        """
        A list of strings containing different command line arguments
        to be used with `gmx mdrun` during benchmarking runs. Every
        entry in the list results in a separate benchmarking run.
        Calling this function before using `write_scripts` is optional.

        Parameters
        ----------
        args_list: List[str]
            A list of strings containing any number of command line arguments to `gmx mdrun`

        """
        assert isinstance(args_list, list) and all([isinstance(v, str) for v in args_list])
        self._run_args = args_list

    def set_mdp_path(self, path: str) -> None:
        """
        Set the path given by the --mdp argument in the calling script. Used
         only to write the add_executable script. Calling this function before using
        `write_scripts` is optional.

        Parameters
        ----------
        path: str
            A string containing the `--mdp` argument

        """
        self._mdp_path = path

    def set_system(
            self, system_name: str, nreps: int, multiply: float, num_system_sizes: int,
            min_boxsize: np.ndarray(shape=(1, 3), dtype=float)) -> None:
        """
        Set details about systems to be ran, namely
          * the name of the system,
          * the number of repetitions for each benchmarking run (to improve statistics),
          * if the system size should be increased: by how much,
          * if the system size should be increased: how many intermediate steps to take,
          * if the system size should be increased: what the minimum box size is.
        Calling this function before using `write_scripts` is mandatory.

        Parameters
        ----------
        system_name: str
            The name of the system, equal to the name of the input files in the respective folder
        nreps: int >= 1
            The number of repetitions for each benchmarking run
        multiply: float > 1 or None
            The maximum box size relative to the minimal box size, if the simulations are to be
            repeated at increasing system size
        num_system_sizes: int >= 2
            If `multiply` is not None, the total number of systems simulated at different box
            sizes, including the minimal and the maximal box
        min_boxsize: np.ndarray(shape=(1,3), dtype=float)
            If `multiply` is not None, the minimal box size

        """
        assert nreps >= 1
        assert multiply is None or (isinstance(multiply, float) and multiply > 1)
        assert multiply is None or (isinstance(num_system_sizes, int) and num_system_sizes >= 2)
        assert multiply is None or isinstance(min_boxsize, np.ndarray)
        self._system_name = system_name
        self._nreps = nreps
        self._multiply = multiply
        self._num_system_sizes = num_system_sizes
        self._min_box = min_boxsize
        self._systems_set = True

    @staticmethod
    def _eq_script_header() -> List:
        """
        (Internal use only)
        Prepares the header of the equilibration file, including help printing

        Returns
        -------
        List
            The header
        """
        return list([
            '#!/bin/bash',
            '',
            'function print_help() {',
            '    echo "USAGE: equilibrate.sh [-h | --help]"',
            '    echo ""',
            '    echo "Run equilibration (needs to be ran only once)."',
            '    echo ""',
            '    echo "ARGUMENTS:"',
            '    echo "  -h, --help            show help message and exit"',
            '    return',
            '}',
            '',
            'for var in "$@"; do',
            '    if [ "$var" == "-h" ] || [ "$var" == "--help" ]; then',
            '        print_help',
            '        exit 0',
            '    else',
            '        echo "ERROR: Unrecognized argument: $var"',
            '        print_help',
            '        exit 1',
            '    fi',
            'done',
            ''])

    def _run_script_header(self) -> List:
        """
        (Internal use only)
        Prepares the header of the run file, including help printing and an argument parser

        Returns
        -------
        List
            The header
        """
        return list([
            '#!/bin/bash',
            '',
            'EXECUTABLES=(all {:s})'.format(' '.join(self._gmx_names)),
            '',
            'function print_help() {',
            '    echo "USAGE: run.sh [-h | --help]"',
            '    local exestring=$(printf " [%s]" "${EXECUTABLES[@]}")',
            '    echo "             $exestring"',
            '    echo ""',
            '    echo "Run benchmarking with one, several or all exectuables."',
            '    echo "Requires equilibration to have ran previously."',
            '    echo ""',
            '    echo "ARGUMENTS:"',
            '    echo "  -h, --help            show help message and exit"',
            '    echo "  all                   run all exectuables"',
            '    for exe in "${EXECUTABLES[@]:1}"; do',
            '        printf "  %-20s  run executable %s\n" $exe $exe',
            '    done',
            '    return',
            '}',
            '',
            'function contains_element() {',
            '    local var match="$1"',
            '    shift',
            '    for var in "$@"; do',
            '        [ "$var" == "$match" ] && return 0',
            '    done',
            '    return 1',
            '}',
            '',
            'if [ $# -le 0 ]; then',
            '    echo "ERROR: No executable chosen."',
            '    print_help',
            '    exit 1',
            'fi',
            '',
            'for var in "$@"; do',
            '    if [ "$var" == "-h" ] || [ "$var" == "--help" ]; then',
            '        print_help',
            '        exit 0',
            '    fi',
            '    contains_element "$var" "${EXECUTABLES[@]}"',
            '    if [ $? -ne 0 ]; then',
            '        echo "ERROR: Unrecognized argument: $var"',
            '        print_help',
            '        exit 1',
            '    fi',
            'done',
            ''])

    @staticmethod
    def _ana_script_header() -> List:
        """
        (Internal use only)
        Prepares the header of the analysis file, including help printing

        Returns
        -------
        List
            The header
        """
        return list([
            '#!/bin/bash',
            '',
            'function print_help() {',
            '    echo "USAGE: analyze.sh [-h | --help]"',
            '    echo ""',
            '    echo "Run analysis for all executables. All executables need to have ran previously."',
            '    echo ""',
            '    echo "ARGUMENTS:"',
            '    echo "  -h, --help            show help message and exit"',
            '    return',
            '}',
            '',
            'for var in "$@"; do',
            '    if [ "$var" == "-h" ] || [ "$var" == "--help" ]; then',
            '        print_help',
            '        exit 0',
            '    else',
            '        echo "ERROR: Unrecognized argument: $var"',
            '        print_help',
            '        exit 1',
            '    fi',
            'done',
            ''])

    @staticmethod
    def _add_script_header() -> List:
        """
        (Internal use only)
        Prepares the header of the add_executable file, including help printing

        Returns
        -------
        List
            The header
        """
        return list([
            '#!/bin/bash',
            '',
            'function print_help() {',
            '    echo "USAGE: add_executable.sh [-h | --help] gmx_exe exe_name"',
            '    echo ""',
            '    echo "Add an executable to the run and the analysis files."',
            '    echo ""',
            '    echo "ARGUMENTS:"',
            '    echo "  gmx_exe               A GROMACS executable to be added to the scripts."',
            '    echo "  exe_name              A unique name for the new executable."',
            '    echo "  -h, --help            show help message and exit"',
            '    return',
            '}',
            '',
            'if [ $# -le 0 ]; then',
            '    echo "ERROR: No executable chosen."',
            '    print_help',
            '    exit 1',
            'fi',
            '',
            'for var in "$@"; do',
            '    if [ "$var" == "-h" ] || [ "$var" == "--help" ]; then',
            '        print_help',
            '        exit 0',
            '    fi',
            'done',
            '',
            'if [ $# -le 1 ]; then',
            '    echo "ERROR: Not enough arguments."',
            '    print_help',
            '    exit 1',
            'fi',
            '',
            'if [ $# -gt 2 ]; then',
            '    echo "ERROR: Too many arguments."',
            '    print_help',
            '    exit 1',
            'fi',
            '',
            'NEWEXE=$1',
            'NEWNAME=$2',
            ''])

    def write_scripts(self) -> None:
        """
        Writes four scripts into the previously defined directory:
          * equilibrate.sh
          * run.sh
          * analyze.sh
          * add_executable.sh

        """
        assert self._directory is not None and self._systems_set
        minimize = True

        do_multiply = self._multiply is not None

        if do_multiply:
            eq_script = self._eq_script_header()
    
            eq_script.append('MAINDIR={:s}\n'.format(self._directory))
            eq_script.append('cd $MAINDIR\n')
    
            eq_script.append('# Using first executable only for equilibration')
            eq_script.append('GMXEXE={:s}'.format(self._gmx_exe[0]))
            eq_script.append('mkdir -p eq')
            eq_script.append('cd eq\n')

            for nsystem, multiplier in enumerate(np.linspace(1, self._multiply, self._num_system_sizes)):
                eq_script.append('## System size {}'.format(nsystem+1))
                eq_script.append('mkdir -p benchmark_{}'.format(nsystem+1))
                eq_script.append('cd benchmark_{}\n'.format(nsystem+1))
                eq_script.append('# System creation & equilibration')
                eq_script.append('cp ../../input/{:s}.top system.top'.format(self._system_name))
                box = self._min_box * multiplier
                eq_script.append('${{GMXEXE}} solvate -box {0:f} {1:f} {2:f} -cs ../../input/{3:s}.gro '
                                 '-p system.top -o start.gro'.format(*box, self._system_name))
                if minimize:
                    eq_script.append('${{GMXEXE}} grompp -f ../../input/{:s}_min.mdp -c start.gro -p system.top '
                                     '-o min.tpr -po min.mdp'.format(self._system_name))
                    eq_script.append('${{GMXEXE}} mdrun -s min.tpr -deffnm min\n'.format())
                    eq_script.append('${{GMXEXE}} grompp -f ../../input/{:s}_eq.mdp -c min.gro -p system.top '
                                     '-o eq.tpr -po eq.mdp'.format(self._system_name))
                else:
                    eq_script.append('${{GMXEXE}} grompp -f ../../input/{:s}_eq.mdp -c start.gro -p system.top '
                                     '-o eq.tpr -po eq.mdp'.format(self._system_name))
                eq_script.append('${{GMXEXE}} mdrun -s eq.tpr -deffnm eq {:s}'.format(self._eq_args))
                eq_script.append('cd ..\n')

            eq_script.append('cd $MAINDIR\n')

            with open(os.path.join(self._directory, 'equilibrate.sh'), 'w') as eq_script_file:
                for line in eq_script:
                    eq_script_file.write(line + '\n')

        run_script = self._run_script_header()
        run_script.append('MAINDIR={:s}\n'.format(self._directory))
        run_script.append('cd $MAINDIR\n')

        for exe, name in zip(self._gmx_exe, self._gmx_names):
            run_script.append('### GMX EXE "{:s}"'.format(name))
            run_script.append('if contains_element "{:s}" "$@" || contains_element "all" "$@"; then'.format(name))
            indent = ' '*4
            run_script.append(indent + 'GMXEXE={:s}'.format(exe))
            run_script.append(indent + 'mkdir -p {:s}'.format(name))
            run_script.append(indent + 'cd {:s}\n'.format(name))

            if not do_multiply:
                run_script.append(indent + 'mkdir -p benchmark_1')
                run_script.append(indent + 'cd benchmark_1\n')
                for nrun, run_args in enumerate(self._run_args):
                    run_script.append(indent + '# Run options {:d}'.format(nrun + 1))
                    if run_args:
                        run_script.append(indent + '# {:s}'.format(run_args))
                    run_script.append(indent + '${{GMXEXE}} grompp '
                                               '-f ../../input/{0:s}_run.mdp -c ../../input/{0:s}.gro '
                                               '-p ../../input/{0:s}.top -o run_{1:d}.tpr -po run_{1:d}.mdp'.format(
                                                   self._system_name, nrun + 1))
                    run_script.append(indent + '${{GMXEXE}} mdrun -s run_{0:d}.tpr -deffnm run_{0:d} '
                                      '{1:s} {2:s}\n'.format(nrun + 1, self._mdrun_args, run_args))

            else:
                for nsystem, multiplier in enumerate(np.linspace(1, self._multiply, self._num_system_sizes)):
                    run_script.append(indent + '## System size {}'.format(nsystem+1))
                    run_script.append(indent + 'mkdir -p benchmark_{}\n'.format(nsystem+1))
                    run_script.append(indent + 'cd benchmark_{}'.format(nsystem+1))
                    run_script.append(indent + '${{GMXEXE}} grompp -f ../../input/{0:s}_run.mdp '
                                      '-c ../../eq/benchmark_{1:d}/eq.gro -p ../../eq/benchmark_{1:d}/system.top '
                                      '-o run.tpr -po run.mdp\n'.format(self._system_name, nsystem+1))

                    for nrun, run_args in enumerate(self._run_args):
                        run_script.append(indent + '# Run options {:d}'.format(nrun+1))
                        if run_args:
                            run_script.append(indent + '# {:s}'.format(run_args))
                        for n in range(self._nreps):
                            run_script.append(indent + '${{GMXEXE}} mdrun -s run.tpr -deffnm run_{0:d}_{1:d} '
                                              '{2:s} {3:s}'.format(nrun+1, n+1, self._mdrun_args, run_args))
                        run_script.append('')

                    run_script.append(indent + 'cd ..\n')

            run_script.append(indent + 'cd $MAINDIR')
            run_script.append('fi\n')

        with open(os.path.join(self._directory, 'run.sh'), 'w') as run_script_file:
            for line in run_script:
                run_script_file.write(line + '\n')
                
        ana_script = self._ana_script_header()
        ana_script.append('MAINDIR={:s}\n'.format(self._directory))
        ana_script.append('cd $MAINDIR\n')

        ana_script.append('gmxbenchmark-analyze \\')
        indent = ' '*4
        dirs = indent
        for name in self._gmx_names:
            dirs += name + ' '
        ana_script.append(dirs + '\\')
        ana_script.append(indent + '--nsystemsizes {:d} \\'.format(self._num_system_sizes))
        for run_args in self._run_args:
            ana_script.append(indent + '--run_args "\\{:s}" \\'.format(run_args))
        ana_script.append(indent + '--nreps {:d}\n'.format(self._nreps))

        ana_script.append('[ $? == 0 ] && echo \'Analysis successful. Open analysis.html to see the results.\'')

        with open(os.path.join(self._directory, 'analyze.sh'), 'w') as ana_script_file:
            for line in ana_script:
                ana_script_file.write(line + '\n')

        add_script = self._add_script_header()
        add_script.append('MAINDIR={:s}\n'.format(self._directory))
        add_script.append('{')
        indent = ' '*4
        add_script.append(indent + 'gmxbenchmark-prepare ${{MAINDIR}} \\')
        indent += ' '*4
        add_script.append(indent + '--gmx {:s} ${{NEWEXE}} \\'.format(' '.join(self._gmx_exe)))
        add_script.append(indent + '--gmxname {:s}  ${{NEWNAME}} \\'.format(' '.join(self._gmx_names)))
        add_script.append(indent + '--system_name {:s} \\'.format(self._system_name))
        add_script.append(indent + '--nreps {:d} \\'.format(self._nreps))
        if do_multiply:
            add_script.append(indent + '--multiply {:g} \\'.format(self._multiply))
            add_script.append(indent + '--nsystemsizes {:d} \\'.format(self._num_system_sizes))
        if self._mdp_path:
            add_script.append(indent + '--mdp {:s} \\'.format(self._mdp_path))
        if len(self._run_args) >= 1 and self._run_args[0] != '':
            for args in self._run_args:
                add_script.append(indent + '--run_args "\\{:s}" \\'.format(args))
        if self._run_args:
            add_script.append(indent + '--mdrun_args "\\{:s}" \\'.format(self._mdrun_args))
        if self._eq_args:
            add_script.append(indent + '--eq_args "\\{:s}" \\'.format(self._eq_args))
        add_script[-1] = add_script[-1][:-2]
        add_script.append(indent + 'exit')
        add_script.append('}')

        with open(os.path.join(self._directory, 'add_executable.sh'), 'w') as add_script_file:
            for line in add_script:
                add_script_file.write(line + '\n')


def make_scripts(args):
    """
    This function copies input files into the simulation folder, and writes scripts
    for initial equilibration, benchmarking simulations, and analysis.

    """
    parser = argparse.ArgumentParser(
        description='Create scripts to run and analyze benchmarks for GROMACS',
        epilog='full example:\n'
               '  python3 gmxbenchmark-prepare ../runs/ --gmx ~/original/bin/gmx ~/modified/bin/gmx \\\n'
               '      --gmxname original modified --multiply 3 --nsystemsizes 9 \\\n'
               '      --run_args "\\-nt 1" --run_args "\\-nt 4 -ntomp 4" --run_args "\\-nt 4 -ntmpi 4" \\\n'
               '      --mdrun_args "\\-pin on" --eq_args "\\-nt 0"\n\n'
               '  This compares two binaries, "original" and "modified", by increasing the box size of the h2o_settle\n'
               '  system up to 3-fold, in 9 steps. The equilibration is ran on all available threads (-nt 0).\n'
               '  The benchmarking simulations are repeated with three sets of arguments, running the simulation on \n'
               '  1 thread (-nt 1), 4 OMP threads (-nt 4 -ntomp 4) threads, and 4 MPI threads (-nt 4 -ntmpi 4),\n'
               '  respectively. All of these use pinning.',
        prog='gmxbenchmark-prepare',
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('directory', type=str, metavar='dir',
                        help='Directory to setup up the simulation files in.')

    parser.add_argument('--gmx', type=str, metavar='exe', nargs='+',
                        help='Give one or more GROMACS executables to run the benchmarking on.\n'
                             'If argument is not used, trying to use \'gmx\' (which needs to be in the PATH).')
    parser.add_argument('--gmxname', type=str, metavar='name', nargs='+',
                        help='Name the gmx executables defined with `--gmx`, used for output only.\n'
                             'If argument is not used, executables will simply be numbered.\n'
                             'If argument is used, length must be equal to number of gmx executables.')
    parser.add_argument('--system_name', type=str, metavar='system', default='h2o_settle',
                        help='The name of the system. Default: h2o_settle')
    parser.add_argument('--nreps', type=int, metavar='N', default=5,
                        help='The number of times the benchmarking runs are repeated, with the timings being\n'
                             'averaged over all repetitions. Default: 5.')
    parser.add_argument('--multiply', type=float, metavar='m', default=None,
                        help='If set, use `gmx solvate` to increase base system size to `m` times\n'
                             'the original box side length. Use `--nsystemsizes` to define\n'
                             'intermediate steps. Default: off.')
    parser.add_argument('--nsystemsizes', type=int, metavar='N', default=2,
                        help='If `--multipy` is set, define how many systems are simulated.\n'
                             'Default: 2 (only original and maximum box size).')

    parser.add_argument('--mdp', type=str, metavar='input.mdp', default=None,
                        help='Give an mdp file allowing to set additional GROMACS options.')
    parser.add_argument('--run_args', type=str, metavar='"\\-gmx_arg X"', action='append',
                        help='Every use of --run_args will prompt an independent benchmarking run of `gmx mdrun`\n'
                             'with the given arguments. Note: Prepend a \'\\\'.\n'
                             'Example: `--run_args "\\-nt 4 -ntomp 4" --run_args "\\-nt 4 -ntmpi 4"`\n'
                             '         will start two runs, one with 4 OMP, the other with 4 MPI threads.')
    parser.add_argument('--eq_args', type=str, metavar='"\\-gmx_arg X"',
                        help='Arguments used for equilibration runs of `gmx mdrun`. Note: Prepend a \'\\\'.')
    parser.add_argument('--mdrun_args', type=str, metavar='"\\-gmx_arg X"',
                        help='Arguments used for every benchmarking run of `gmx mdrun`. Note: Prepend a \'\\\'.')

    args = parser.parse_args(args)
    
    if args.run_args:
        for idx in range(len(args.run_args)):
            args.run_args[idx] = args.run_args[idx][1:]
    if args.eq_args:
        args.eq_args = args.eq_args[1:]
    if args.mdrun_args:
        args.mdrun_args = args.mdrun_args[1:]

    no_multiplication = not args.multiply

    # prepare directory
    input_directory = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'input')
    simulation_directory = os.path.abspath(args.directory)
    pathlib.Path(simulation_directory).mkdir(parents=True, exist_ok=True)
    simulation_input_directory = os.path.join(simulation_directory, 'input')
    pathlib.Path(simulation_input_directory).mkdir(exist_ok=True)

    # prepare mdp files
    shutil.copy2(os.path.join(input_directory, args.system_name + '_min.mdp'), simulation_input_directory)
    shutil.copy2(os.path.join(input_directory, args.system_name + '_eq.mdp'), simulation_input_directory)
    shutil.copy2(os.path.join(input_directory, args.system_name + '_run.mdp'), simulation_input_directory)
    if args.mdp is not None:
        mdp_eq_file = os.path.join(simulation_input_directory, args.system_name + '_eq.mdp')
        mdp_run_file = os.path.join(simulation_input_directory, args.system_name + '_run.mdp')
        mdp_dict_eq = read_mdp(mdp_eq_file)
        mdp_dict_run = read_mdp(mdp_run_file)
        mdp_dict_eq = read_mdp(args.mdp, mdp_dict_eq)
        mdp_dict_run = read_mdp(args.mdp, mdp_dict_run)
        write_mdp(mdp_eq_file, mdp_dict_eq)
        write_mdp(mdp_run_file, mdp_dict_run)

    # copy topology
    if no_multiplication:
        shutil.copy2(os.path.join(input_directory, args.system_name + '.top'), simulation_input_directory)
    else:
        strip_topology(
            os.path.join(input_directory, args.system_name + '.top'),
            os.path.join(simulation_input_directory, args.system_name + '.top')
        )

    # copy coordinate files
    shutil.copy2(os.path.join(input_directory, args.system_name + '.gro'), simulation_input_directory)
    min_boxsize = get_box_size(os.path.join(input_directory, args.system_name + '.gro'))

    script = Script()

    script.set_directory(simulation_directory)
    script.set_system(args.system_name, args.nreps, args.multiply, args.nsystemsizes, min_boxsize)

    if args.gmx:
        gmx_exe = args.gmx
    else:
        gmx_exe = ['gmx']

    gmx_names = []
    if args.gmxname:
        assert len(gmx_exe) == len(args.gmxname)
        assert all(('/' not in name for name in args.gmxname))  # assert no paths were given as names
        gmx_names = args.gmxname
    else:
        for n in enumerate(gmx_exe):
            gmx_names.append('GMX_{:d}'.format(n))

    if len(gmx_names) != len(set(gmx_names)):  # assert no duplicate names
        raise RuntimeError('Names (--gmx_names) must be unique.')
    script.set_gmx(gmx_exe, gmx_names)

    if args.run_args:
        script.set_run_args(args.run_args)
    if args.eq_args:
        script.set_eq_args(args.eq_args)
    if args.mdrun_args:
        script.set_mdrun_args(args.mdrun_args)

    if args.mdp is not None:
        script.set_mdp_path(args.mdp)

    script.write_scripts()

    return 0


def main():
    sys.exit(make_scripts(sys.argv[1:]))


if __name__ == '__main__':
    main()
