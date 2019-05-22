**********************
GROMACS benchmark tool
**********************

This tool aims at evaluating performance of the GROMACS molecular
simulation engine (www.gromacs.org). In the short term, this is mostly
thought to be a help for code development, in assessing how changes to
the code influence its performance. In the longer term, this tool could
evolve to also allow the evaluation of hardware and parallelization
schemes based on more complex, possibly user-submitted systems.

Contents
========

* `Installation`_
* `Usage`_

  * `Prepare scripts`_
  * `Run scripts`_
  * `Systems`_

* `Analysis`_


Installation
============

Development version
-------------------

The latest version is available on `github`_. It can be installed directly using ``pip``
::

   pip install git+https://github.com/ptmerz/gmxbenchmark.git

Alternatively, the source code can also be downloaded from `github`_ and then installed
using the provided setup file,
::

   python3 setup.py install


.. _`github`: https://github.com/ptmerz/gmxbenchmark

Usage
=====

Prepare scripts
---------------

``gmxbenchmark-prepare`` allows users to prepare equilibration, run, and analysis scripts according to
their needs. It is ran directly from the command line.

Command-line help
~~~~~~~~~~~~~~~~~
::

   $ gmxbenchmark-prepare -h
   usage: gmxbenchmark-prepare [-h] [--gmx exe [exe ...]]
                               [--gmxname name [name ...]] [--system_name system]
                               [--nreps N] [--multiply m] [--nsystemsizes N]
                               [--mdp input.mdp] [--run_args "\--gmx_arg X"]
                               [--eq_args "\--gmx_arg X"]
                               [--mdrun_args "\--gmx_arg X"]
                               dir

   Create scripts to run and analyze benchmarks for GROMACS

   positional arguments:
     dir                   Directory to setup up the simulation files in.

   optional arguments:
     -h, --help            show this help message and exit
     --gmx exe [exe ...]   Give one or more GROMACS executables to run the benchmarking on.
                           If argument is not used, trying to use 'gmx' (which needs to be in the PATH).
     --gmxname name [name ...]
                           Name the gmx executables defined with `--gmx`, used for output only.
                           If argument is not used, executables will simply be numbered.
                           If argument is used, length must be equal to number of gmx executables.
     --system_name system  The name of the system. Default: h2o_settle
     --nreps N             The number of times the benchmarking runs are repeated, with the timings being
                           averaged over all repetitions. Default: 5.
     --multiply m          If set, use `gmx solvate` to increase base system size to `m` times
                           the original box side length. Use `--nsystemsizes` to define
                           intermediate steps. Default: off.
     --nsystemsizes N      If `--multipy` is set, define how many systems are simulated.
                           Default: 2 (only original and maximum box size).
     --mdp input.mdp       Give an mdp file allowing to set additional GROMACS options.
     --run_args "\-gmx_arg X"
                           Every use of --run_args will prompt an independent benchmarking run of `gmx mdrun`
                           with the given arguments. Note: Prepend a '\'.
                           Example: `--run_args "\-nt 4 -ntomp 4" --run_args "\-nt 4 -ntmpi 4"`
                                    will start two runs, one with 4 OMP, the other with 4 MPI threads.
     --eq_args "\-gmx_arg X"
                           Arguments used for equilibration runs of `gmx mdrun`. Note: Prepend a '\'.
     --mdrun_args "\-gmx_arg X"
                           Arguments used for every benchmarking run of `gmx mdrun`. Note: Prepend a '\'.

   full example:
     python3 gmxbenchmark-prepare ../runs/ --gmx ~/original/bin/gmx ~/modified/bin/gmx \
         --gmxname original modified --multiply 3 --nsystemsizes 9 \
         --run_args "\-nt 1" --run_args "\-nt 4 -ntomp 4" --run_args "\-nt 4 -ntmpi 4" \
         --mdrun_args "\-pin on" --eq_args "\-nt 0"

     This compares two binaries, "original" and "modified", by increasing the box size of the h2o_settle
     system up to 3-fold, in 9 steps. The equilibration is ran on all available threads (-nt 0).
     The benchmarking simulations are repeated with three sets of arguments, running the simulation on
     1 thread (-nt 1), 4 OMP threads (-nt 4 -ntomp 4) threads, and 4 MPI threads (-nt 4 -ntmpi 4),
     respectively. All of these use pinning.

Details
~~~~~~~

**Run directory, executables and system name:** The only positional argument is the run directory, a
directory to which run input files (parameter, topology, starting structures) and run scripts are written
to. Typical usage of the benchmarking suite will consist of performance comparison of two or more
GROMACS executables. The executables used in this comparison can be specified using the ``--gmx`` flag.
For reporting purposes, the executables are named. The (unique) names can be chosen with the ``--gmxname``
argument. If no names are chosen, the executables are named ``EXE_1``, ``EXE_2``, etc. The system to be
simulated is chosen using the ``--system_name`` argument. See `Systems`_ for a list of available systems.

A simple call to ``gmxbenchmark-prepare`` to set up simulations comparing the original GROMACS binary to a
modified version using a simple water system (using settle to constrain bonds) could hence look like
::

   gmxbenchmark-prepare my/simuldir \
       --gmx ~/original/bin/gmx ~/modified/bin/gmx \
       --gmxname original modified \
       --system_name h2o_settle

**Repetitions:** The ``--nreps`` argument allows to set how often the benchmarking is repeated. Repeating
the benchmarking is encouraged to produce more reliable results. If ``--nreps`` is not explicitly set, the
benchmarking is repeated five times.

**Additional mdp options:** To test specific setups, it might be useful to define (or overwrite) specific
parameters of the GROMACS input file. Options can be given via a ``mdp`` file and the ``--mdp`` argument.
Parameters found in this file are added to the standard settings defined in the input files. If options
are found in both files, the options found via user-input are prioritized.

**Multiple system sizes:** To test performance at different system sizes in a systematic way,
``gmxbenchmark-prepare`` allows to defined the ``--multiply m`` argument. This takes the box size of the
chosen system, and creates a system whose box edges are ``m`` times the original size. ``--nsystemsizs``
allows to define the number of systems created. The standard is 2 - the original system size and the
multiplied box. Choosing a higher number creates intermediate system sizes. The boxes generated in that
way are filled using ``gmx solvate``. This has only been tested with systems consisting of a single
type of molecules. Note that creating systems in this way requires an equilibration step (see
`Equilibration`_ below), while non-mulitplied systems are assumed to be equilibrated and ready for
benchmarking.

Example:
::

   gmxbenchmark-prepare my/simuldir \
       --gmx ~/original/bin/gmx ~/modified/bin/gmx \
       --gmxname original modified \
       --system_name h2o_settle \
       --multiply 3 --nsystemsizes 9

This would use the original ``h2o_settle`` system (whose box size is 2.2nm x 2.2nm x 2.2nm) to create
nine systems of box sizes ranging from (2.2nm)^3 to (6.6nm)^3. The box lengths of the systems are
equidistant, i.e. 2.2nm, 2.75nm, 3.3nm, 3.85nm, 4.4nm, 4.95nm, 5.5nm, 5.05nm, 6.6nm. The number of
atoms in the system would range from 1044 (original size) to 28188.

**Run arguments:** There are three types of arguments which can be added to the ``gmx mdrun`` commands
used. ``--eq_args`` are used for equilibration runs only. ``--mdrun_args`` are used for all ``gmx mdrun``
calls used during benchmarking. Most importantly, ``--run_args``, which can be used repeatedly, define
different benchmarking settings. This allows to compare the executables in a range of different
setups. For example,
::

   gmxbenchmark-prepare my/simuldir \
       --gmx ~/original/bin/gmx ~/modified/bin/gmx \
       --gmxname original modified \
       --system_name h2o_settle \
       --multiply 3 --nsystemsizes 9 \
       --run_args "\-nt 1" --run_args "\-nt 4 -ntomp 4" --run_args "\-nt 4 -ntmpi 4 -dlb yes" \
       --mdrun_args "\-pin on" --eq_args "\-nt 0"

will repeat the benchmarking simulations defined earlier with three different sets of ``gmx mdrun``,
running on a single thread, on 4 OpenMP threads, and on 4 MPI threads. All of them will use pinning.
The equilibrations, on the other hand, will run on all avaiable resources.

.. note:: Having repetitions of the same benchmarking (same executable, same run arguments, same system
   size) differ in whether dynamic load balancing is performed or not is currently not handled by the
   analysis script. It is therefore advised to explicitly turn load balancing on or off (when applicable).

Run scripts
-----------
``gmxbenchmark-prepare`` writes four scripts to the selected directory, ``equilibrate.sh``, ``run.sh``,
``add_executable.sh``, and ``analyze.sh``. These scripts are described below.

Equilibration
~~~~~~~~~~~~~
Equilibration is only performed if the ``--multiply`` option was chosen, as other systems are expected
to be equilibrated. ``equilibrate.sh`` doesn't take any arguments, and will first create systems of
different sizes (using ``gmx solvate``) and minimize them, before running a longer equilibration run.
The results of these equilibration runs are found in the ``eq/`` folder. The equilibration is only ran
with one executable (the first one given in the ``--gmx`` arguments). It does hence not need to be reran
if new executables are added (see `Add executables`_ below), or if the benchmarking is re-ran after code
changes.

Benchmarking run
~~~~~~~~~~~~~~~~
To benchmark the executables, the systems chosen are ran repeatedly with all executables to gather
data for later analysis. ``run.sh`` allows to specify which executable to run (``all`` or one or more
executables by their names defined through ``--gmxname``). Running executables separately allows,
for example, to re-run only a specific executable after recompiling it, or run only a newly added
executable (see `Add executables`_ below). Use
::

   run.sh -h

to see all options.

Add executables
~~~~~~~~~~~~~~~
To add an executable later on, one can run ``gmxbenchmark-prepare`` again with the same arguments as the
first time, only adding the additional executable (and, if names were chosen, an additional name).
``add_executable.sh`` simplifies this process by only taking an executable and a name as input and rerunning
``gmxbenchmark-prepare`` with the previously chosen arguments. Note that since the equilibration, if needed,
is only ran with one executable, it doesn't need to be reran upon addition of an executable. Also note that
``run.sh`` allows to run specific executables only, so it also isn't necessary to rerun benchmarking for all
executables.

Analyze
~~~~~~~
``gmxbenchmark-analyze`` reads the timing reported in the log files of runs prepared by ``gmxbenchmark-prepare``,
and reports the mean and standard deviations in an html file named ``analysis.html``. ``gmxbenchmark-analyze``
requires some arguments mirroring the ones used for ``gmxbenchmark-prepare`` to be able to retrieve the
results. ``analyze.sh`` simplifies this by providing a simple script with no input arguments needed.

Systems
-------
The following systems are currently available:

``h2o_settle``
~~~~~~~~~~~~~~
*(This is the default system)*

Pure water box, 348 molecules, 1044 atoms, box size 2.2nm x 2.2nm x 2.2nm. The topology prescribes constraints
using settle. Uses reaction field with ``epsilon-rf = 0``. The equilibration script uses thermostatting, while
the benchmarking simulation is ran in NVE. Note: Changing the ensemble sampled can easily be done using the
``--mdp`` option of ``gmxbenchmark-prepare``.

``h2o``
~~~~~~~
Pure water box, 348 molecules, 1044 atoms, box size 2.2nm x 2.2nm x 2.2nm. The topology prescribes no
constraints. Uses reaction field with ``epsilon-rf = 0``. The equilibration script uses thermostatting, while
the benchmarking simulation is ran in NVE. Note: Changing the ensemble sampled can easily be done using the
``--mdp`` option of ``gmxbenchmark-prepare``.

Contributing new systems
~~~~~~~~~~~~~~~~~~~~~~~~
New systems should come with an equilibrated starting structure, a topology file, a ``mdp``-parameter
file for the benchmarking run, and a short description. If the system is suitable to be ran at different
sizes using the ``--multiply`` option, ``mdp``-parameter files for minimization and equilibration should
be added.

Analysis
========
Running the analysis script will create ``analysis.html`` which can be opened using any browser. An example
of such a result file con be found
`here <https://htmlpreview.github.io/?https://github.com/ptmerz/gmxbenchmark/blob/master/examples/analysis.html>`_.
The file consists of a few different sections, namely

* **Key:** Describes the executable, the benchmarking run arguments, and the system sizes used, and assigns
  each a key used in the remainder of the file. In case of the example, two builds of the same GROMACS version
  are compared. They were compared in three different run scenarios, namely on a single thread (*arg1*), on 4
  OpenMP threads (*arg2*), and on 4 MPI threads (*arg3*), and using 9 different system sizes.
* **Plot:** Next, the core time normalized by the total number of atoms is plotted. A click on the plot allows
  to view a larger version of the figure.
* **Core times:** The data plotted in the figure above are given in numbers. A click on the blue ⓘ icon leads
  to detailed timings for the specific simulation (see below for details).
* **Compare timings:** This allows to easily compare the timings (again, see below for details) of any choice
  of argument and system size across different executables by displaying the timing blocks next to each other.
* **Timings:** The timings, finally, consist of the mean values and standard deviations over the benchmarking
  repetitions of the timing blocks found in the GROMACS log files. They allow to investigate further where
  differences between executables arise.
