from typing import Dict
import re
import os
from collections import defaultdict
import pandas as pd


class BrokenLogError(Exception):
    pass


class LammpsLog():
    filepath: str
    product_idx: int
    """
    Parse a Lammps log file.
    """
    # floating-point number with scientific notation:
    # numbers such as .12 are not expected.
    _real_num = re.compile(r"[-+]?\d+(?:\.\d+)?(?:[eE][\-\+]?\d+)?")

    def __init__(self, filepath: str, product_idx: int) -> None:
        self.filepath = filepath
        self.filename, _ = os.path.splitext(self.filepath)
        self.filename: str = self.filename.split("/")[-1]
        self.product_idx = product_idx
        with open(self.filepath, 'r') as logfile:
            self.log_txt = logfile.read()
        self._find_chunks()
        self._thermo_kws()

    def __str__(self) -> str:
        lammps_log: str = f"""
        LAMMPS log file: '{self.filename}'
        """
        return lammps_log

    def __repr__(self) -> str:
        return (f"LAMMPS log file('{self.filename}'")

    def _find_chunks(self):
        """
        Find the positions of the head(s) and foot(s) of the thermo chunk(s),
        and the position of the wall time info in `self.log_txt`.

        In a given Lammps log find, depending on the numbers calls of Lammps
        'run' command, there is one or more chunk(s) of the thermodynamic
        outputs, mpi statistics, minimazation report, and the like. Each
        thermodynamicchunk starts with a head line about MPI memory usage and
        ends with a foot line about the run time. Other information than the
        thermodynamic outputs come before the head of the first chunk, between
        the head and foot of two consecutive chunks, or between the foot of
        the last chunk and the 'Total wall time' line.

        Raises
        ------
        BrokenLogError
            The numbers of headers and footers are not equal.
        """
        # See these references for the patterns used to find the chunks:
        # https://github.com/lammps/lammps/blob/2b8d6fc4d93f6c7dcce870ed134e6d3ac45f47b7/src/finish.cpp#L129
        # https://github.com/lammps/lammps/blob/27b37ebad04511e49d9e3c0303390bf811194095/src/output.cpp#L1027
        # https://github.com/lammps/lammps/blob/dd8c1df9c24c623f910801dee6e4b34d9d8d2065/src/lammps.cpp#L757
        head_pat = re.compile(
            "(Memory usage per processor|Per MPI rank memory)(.*?)\n"
            )
        self.thermo_heads = list(re.finditer(head_pat, self.log_txt))
        foot_pat = re.compile("Loop time of(.*?)\n")
        self.thermo_foots = list(re.finditer(foot_pat, self.log_txt))
        wall_time_pat = re.compile("Total wall time(.*?)\n")
        self.wall_time_txt = list(re.finditer(wall_time_pat, self.log_txt))
        if len(self.thermo_heads) != len(self.thermo_foots):
            raise BrokenLogError(
                "The number of headers and footers are not equal."
            )
        if self.wall_time_txt == []:
            raise BrokenLogError(
                "No 'Total wall time' found."
            )

    def _thermo_kws(self):
        """
        Parse the first thermo chunk to extract the thermodynamic keywords
        from the first few lines of the first chunk.

        LAMMPS has foor different thermo styles: multi, one, yaml, and custom.
        While the header can be directly defined for multi, one, and yaml
        styles, a general approach is used based on regex to handle all thermo
        styles including custom.

        In the one, yaml, and custom styles, the first line of the chunk is
        the thermodynamic keywords and the rest of lines are the values dumped
        in every few steps;hereafter, each of these line is called a section.
        However, each section in the multi style is composed of 5 lines. In
        contrary to other styles, a section in multi style has both keywords
        and the values, so we must a different section splitter than '\n' for
        the mutli style.

        References
        ----------
        Lammps "thermo" source file:
            https://github.com/lammps/lammps/blob/develop/src/thermo.cpp
        """
        # The thermodynamic keywrods are a combination of alphabetic
        # characters and underscore:
        kw_pat = re.compile("\b*([A-Za-z_]{2,})\b*")
        first = self.thermo_heads[self.product_idx].end()
        last = self.thermo_foots[self.product_idx].start()
        chunk = self.log_txt[first:last-1]
        # The thermo keywords in all the thermo styles except the multi style
        # are in the first line of each chunk. In the multi style, the
        # keywords are spreaded in several line. According to the LAMMPS log
        # file, each thermo output has 5 lines in the multi style, so we use
        # the frist 5 lines of the first chunk in the production phase to find
        # thermodynamics. To make approach general, we do in a messy way.
        # There is no themodynamic keywords with a single-character name, so
        # use {2,} to ignore 'e' or 'E' in numbers with the scientific format
        head_lines = ' '.join(chunk.split("\n")[:5])
        words = re.findall(kw_pat, head_lines)
        try:
            words.remove('sec')  # remove 'sec' word in 'multi' style
            self.sec_split = "\n--"
            self.sec_start = 0
        except ValueError:
            self.sec_split = "\n"
            self.sec_start = 1
        self.keywords = words

    def extract_thermo(self):
        """
        Parse thermo chunks to extract the thermodynamic variable and create
        a Pandas dataframe from the extracted thermodynamic keywords and data.

        Lammps has foor different thermo styles: multi, one, yaml, and custom.
        While the header can be directly defined for multi, one, and yaml
        styles, a general approach is used based on regex to handle all thermo
        styles including custom.

        References
        ----------
        Lammps "thermo" source file:
            https://github.com/lammps/lammps/blob/develop/src/thermo.cpp
        """
        thermo_data = []
        for head, foot in zip(
            self.thermo_heads[self.product_idx:],
            self.thermo_foots[self.product_idx:]
        ):  # loop over thermo chunks
            sections = \
                self.log_txt[head.end():foot.start()-1].split(self.sec_split)
            for section in sections[self.sec_start:]:
                values = re.findall(self._real_num, section)
                thermo_data.append(list(map(float, values)))
        self.thermo = pd.DataFrame(thermo_data, columns=self.keywords)

    def extract_run_stat(self):
        """
        Parse a LAMMPS log file to extract the performance information about a
        given run (chunk) in a full run of a simulation.

        References
        ----------
        https://github.com/lammps/lammps/blob/2b8d6fc4d93f6c7dcce870ed134e6d3ac45f47b7/src/finish.cpp#L129
        """
        run_data: Dict[str, list] = defaultdict(list)
        wall_time = {}
        run_foots = self.thermo_heads[1:] + self.wall_time_txt
        for idx, (head, foot) in enumerate(zip(self.thermo_foots, run_foots)):
            run_data['loop_idx'].append(idx+1)  # loop/chunk idx
            lines = self.log_txt[head.start():foot.end()].split("\n")
            iline = 0
            n_cores = 0
            n_atoms = 0
            while iline < len(lines):
                line = lines[iline]
                if line.startswith('Loop time'):
                    values = re.findall(self._real_num, line)
                    run_data['loop_time'].append(float(values[0]))
                    n_cores = int(values[1])
                    run_data['n_cores'].append(n_cores)
                    run_data['loop_timesteps'].append(int(values[2]))
                    n_atoms = int(values[3])
                    run_data['n_atoms'].append(n_atoms)
                if line.startswith('Performance:'):
                    # handling situation where "atom-step/s" is missed:
                    values = line.split("timesteps/s")[0]
                    values = re.findall(self._real_num, values)
                    if len(values) == 2:  # lj unit
                        run_data['tau_day'].append(float(values[0]))
                        run_data['timestep_sec'].append(float(values[1]))
                    else:  # other units
                        run_data['tau_day'].append(float(values[0]))
                        run_data['timestep_sec'].append(float(values[2]))
                if re.match(r"^\d{1,2}.?\d+\% CPU", line):
                    values = re.findall(self._real_num, line)
                    run_data['cpu_use'].append(float(values[0]))
                    run_data['mpi_task'].append(int(values[1]))
                    try:
                        run_data['openmp_threads'].append(int(values[2]))
                    except IndexError:  # handing no openmp threads
                        run_data['openmp_threads'].append(0)
                if line.startswith("MPI task timing breakdown:"):
                    # only "avg_time" and "%total" are used:
                    iline = iline + 3
                    line = lines[iline]
                    while len(line) != 0:
                        mpi_kw = line.split("|")[0].strip()
                        values = re.findall(self._real_num, line)
                        values = list(map(float, values))
                        run_data[mpi_kw+'_avg_time'].append(values[0])
                        run_data[mpi_kw+'_total_pct'].append(values[1])
                        iline += 1
                        line = lines[iline]
                if line.startswith('Dangerous builds'):
                    words = line.split()
                    run_data['dangerous_builds'].append(int(words[-1]))
                if line.startswith('Total wall time:'):
                    wall_time['n_cores'] = [n_cores]
                    wall_time['n_atoms'] = [n_atoms]
                    h, m, s = line.split()[-1].strip().split(":")
                    wall_t_hr = (int(h) * 3600 + int(m) * 60 + int(s)) / 3600
                    wall_time['wall_time_hr'] = [round(wall_t_hr, 2)]
                iline += 1
        self.run_stat = pd.DataFrame.from_dict(run_data)
        self.wall_time = pd.DataFrame.from_dict(wall_time)
