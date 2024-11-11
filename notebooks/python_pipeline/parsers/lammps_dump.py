from typing import IO, Tuple, Dict, List
from collections import OrderedDict
import numpy as np
from polyphys.utilizer import openany_context, InputT

class Snapshot:
    """
    Read a single snapshot from a dump `file` object.
    """
    def __init__(self, file: InputT) -> None:
        self.file = file
        # time, n_atoms, and n_props are default to invalid numbers.
        self.time = -1  # time stamp
        self.n_atoms = 0  # number of atoms
        self.boxstr = ''  # box information
        self.is_scaled = -1  # unknown status for coordinates
        self.triclinic: bool = False  # default: orthogonal box
        self.xlo = 0.0  # box information
        self.xhi = 0.0
        self.xy = 0.0
        self.ylo = 0.0
        self.yhi = 0.0
        self.xz = 0.0
        self.zlo = 0.0
        self.zhi = 0.0
        self.yz = 0.0
        # names and col indices of per-atom properties:
        self.props: Dict[str, int] = {}
        self.n_props = 0  # number of per atom properties in the dump file
        self.atoms = np.empty((self.n_atoms, self.n_props), dtype=np.float64)
        self._read()

    def _read(self) -> None:
        """
        Read a single snapshot and assigns column names (file must be
        self-describing). Moreover, it changes `self.is_scaled` from -1
        (unknown) to 0 (unscaled) or 1 (scaled), and converts 'xs' and 'xu' to
        'x' in `self.props`.
        """
        _ = self.file.readline()
        # just grab 1st field
        self.time = int(self.file.readline().split()[0])
        _ = self.file.readline()
        self.n_atoms = int(self.file.readline())
        self.boxstr = self.file.readline()
        words = self.boxstr.split()
        self.triclinic = words == 9  # ITEM: BOX BOUNDS
        if self.triclinic:
            self.xlo, self.xhi, self.xy = \
                map(float, self.file.readline().split())
            self.ylo, self.yhi, self.xz = \
                map(float, self.file.readline().split())
            self.zlo, self.zhi, self.yz = \
                map(float, self.file.readline().split())
        else:
            self.xlo, self.xhi = map(float, self.file.readline().split())
            self.ylo, self.yhi = map(float, self.file.readline().split())
            self.zlo, self.zhi = map(float, self.file.readline().split())

        x_flag = y_flag = z_flag = -1
        self.props = {
            c: i for i, c in enumerate(str(self.file.readline().split()[2:]))
            }  # dump per atom props
        props = list(self.props.keys())
        if ("x" in props) or ("xu" in props):
            x_flag = 0
        elif ("xs" in props) or ("xsu" in props):
            x_flag = 1
        elif ("y" in props) or ("yu" in props):
            y_flag = 0
        elif ("ys" in props) or ("ysu" in props):
            y_flag = 1
        elif ("z" in props) or ("zu" in props):
            z_flag = 0
        elif ("zs" in props) or ("zsu" in props):
            z_flag = 1
        if x_flag == 0 and y_flag == 0 and z_flag == 0:
            self.is_scaled = 0
        if x_flag == 1 and y_flag == 1 and z_flag == 1:
            self.is_scaled = 1
        words = self.file.readline().split()
        self.n_props = len(words)
        self.atoms = np.zeros((self.n_atoms, self.n_props), dtype=np.float64)
        self.atoms[0, :] = np.array(list(map(float, words)))  # read first atom
        if self.n_props != len(self.props):
            raise ValueError(
                f"Number of per atom columns '{self.n_props}' is not equal "
                f"to the number of column names '{len(self.props)}'."
                )
        for i in range(1, self.n_atoms):
            self.atoms[i, :] = \
                np.array(list(map(float, self.file.readline().split())))


class Dump:
    """
    Read, write, and manipulate dump files and particle attributes.

    Examples
    --------
    d = dump("dump.one")              read in a single dump.
    d = dump(["dump.1", "dump.2"])	  can be more than one dump.
    d = dump(["dump.1", "dump.2.gz"]) can be gzipped

    `Dump` automatically delete incomplete and duplicate snapshots and
    unscaled coordinates if they are stored in files as scaled.
    """
    def __init__(self, filepath: str) -> None:
        self.filepath = filepath
        self.names: Dict[str, int] = {}
        self.is_scaled = -1  # -1/0/1 mean unknown/unscale/scale coordinates.
        self.snaps: List[Snapshot] = []
        self.n_snaps: int = 0  # total number of snapshots
        self.increment = 1
        self.eof = 0

    def read_all(self) -> None:
        """
        Read all snapshots from each file; test for gzipped files.
        """
        with openany_context(self.filepath) as dfile:
            snap = Snapshot(dfile)
            if self.names == {}:
                self.names = snap.props
                print("Assigned columns: " + self.names2str())
                self.is_scaled = snap.is_scaled
            while snap:
                self.snaps.append(snap)
                snap = Snapshot(dfile)
        # sort entries by timestep, cull duplicates
        self.snaps.sort(key=self.compare_time)
        self.cull()
        self.n_snaps = len(self.snaps)
        print(f"{self.n_snaps} snapshots were read.")
        # select all timesteps and atoms
        # if snapshots are scaled, unscale them
        if ("x" not in self.names.keys()) or \
            ("y" not in self.names.keys()) or \
                ("z" not in self.names.keys()):
            print("Dump scaling status is unknown.")
        elif self.n_snaps > 0:
            if self.is_scaled == 1:
                self.unscale()
            elif self.is_scaled == 0:
                print("dump is already unscaled")
        else:
            print("Dump scaling status is unknown")

    def names2str(self) -> str:
        """
        Convert column names to a string, by column order.

        Returns
        -------
        col_str: str
            The string of column names. `col_str` is repeated in the header
            section of a snapshot, determining per atom properties in that
            snapshot.
        """
        names_res: Dict[int, str] = \
            dict(((v, k) for k, v in self.names.items()))
        names_res = OrderedDict(sorted(names_res.items()))
        col_str = " ".join(names_res.values())
        return col_str

    def cull(self) -> None:
        """
        Delete successive snapshots with duplicate time stamp.
        """
        i = 1
        while i < len(self.snaps):
            if self.snaps[i].time == self.snaps[i-1].time:
                del self.snaps[i]
            else:
                i += 1

    def find_snapshot(self, ts: int) -> int:
        """
        Find a snapshot with timestep `ts`.

        Parameters
        ----------
        ts: int
            The snapshot to be found.

        Returns
        -------
        i: int
            The index of the snapshot in `self.n_snaps` for which timestep is
            equal to `ts`.

        Raises
        ------
        ValueError
            raise error  if `ts` is not found.
        """
        for i in range(self.n_snaps):
            if self.snaps[i].time == ts:
                return i
        raise ValueError(f"No step '{ts}' exists.")

    @staticmethod
    def compare_time(snap: Snapshot) -> int:
        """
        Sort snapshots on time stamp.
        """
        return snap.time

    @staticmethod
    def unit_cell(snap: Snapshot) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute elements of the h-matrix for a tricilinic or orthogonal unit
        cell and its inverse matrix
        Parameters
        ----------
        snap : Snapshot
            Snapshot being unscaled.

        Return
        ------
        h_mat : numpy array
            An array of length 6 containing the elements of h_matrix.
        h_inv : numpy array
            An array of length 6 containing the elements of the inverse of
            the h_matrix.
        """
        xlo_bound = snap.xlo
        xhi_bound = snap.xhi
        ylo_bound = snap.ylo
        yhi_bound = snap.yhi
        zlo_bound = snap.zlo
        zhi_bound = snap.zhi
        xy = snap.xy
        xz = snap.xz
        yz = snap.yz
        xlo = xlo_bound - min((0.0, xy, xz, xy+xz))
        xhi = xhi_bound - max((0.0, xy, xz, xy+xz))
        ylo = ylo_bound - min((0.0, yz))
        yhi = yhi_bound - max((0.0, yz))
        zlo = zlo_bound
        zhi = zhi_bound
        h_mat = np.zeros(6, dtype=np.float64)
        h_mat[0] = xhi - xlo
        h_mat[1] = yhi - ylo
        h_mat[2] = zhi - zlo
        h_mat[3] = yz
        h_mat[4] = xz
        h_mat[5] = xy
        h_inv = np.zeros(6, dtype=np.float64)
        h_inv[0] = 1.0 / h_mat[0]
        h_inv[1] = 1.0 / h_mat[1]
        h_inv[2] = 1.0 / h_mat[2]
        h_inv[3] = yz / (h_mat[1]*h_mat[2])
        h_inv[4] = (h_mat[3]*h_mat[5] - h_mat[1]*h_mat[4]) / \
            (h_mat[0]*h_mat[1]*h_mat[2])
        h_inv[5] = xy / (h_mat[0]*h_mat[1])
        return h_mat, h_inv

    def scale(self, timestep: int = -1) -> None:
        """
        Scale coordinates from box size to [0,1] for all snapshots or for
        snapshot `timestep`, using the h-matrix to treat orthogonal or
        triclinic boxes.

        Parameters
        ----------
        timestep : int, optional
            The timestep of the snapshot being unscaled.
        """
        if timestep == -1:
            print("Scaling dump ...")
            x = self.names["x"]
            y = self.names["y"]
            z = self.names["z"]
            for snap in self.snaps:
                self.scale_one(snap, x, y, z)
        else:
            i = self.find_snapshot(timestep)
            x = self.names["x"]
            y = self.names["y"]
            z = self.names["z"]
            self.scale_one(self.snaps[i], x, y, z)

    def scale_one(self, snap: Snapshot, x: int, y: int, z: int) -> None:
        """
        Scale a single snapshot.

        Parameters
        ----------
        snap : Snapshot
            Snapshot being unscaled.
        x : int
            Column id of the x coordinate.
        y : int
            Column id of the y coordinate.
        z : int
            Column id of the z coordinate.
        """
        # scale coordinates in a single `snap`, using
        if not snap.triclinic:
            xprdinv = 1.0 / (snap.xhi - snap.xlo)
            yprdinv = 1.0 / (snap.yhi - snap.ylo)
            zprdinv = 1.0 / (snap.zhi - snap.zlo)
            atoms = snap.atoms
            if atoms is not None:
                atoms[:, x] = (atoms[:, x] - snap.xlo) * xprdinv
                atoms[:, y] = (atoms[:, y] - snap.ylo) * yprdinv
                atoms[:, z] = (atoms[:, z] - snap.zlo) * zprdinv
        else:
            # Creating h-matrix:
            _, h_inv = self.unit_cell(snap)
            atoms = snap.atoms
            if atoms is not None:
                atoms[:, x] = (atoms[:, x] - snap.xlo)*h_inv[0] + \
                    (atoms[:, y] - snap.ylo)*h_inv[5] + \
                    (atoms[:, z] - snap.zlo)*h_inv[4]
                atoms[:, y] = (atoms[:, y] - snap.ylo)*h_inv[1] + \
                    (atoms[:, z] - snap.zlo)*h_inv[3]
                atoms[:, z] = (atoms[:, z] - snap.zlo)*h_inv[2]

    def unscale(self, timestep: int = -1) -> None:
        """
        Unscale coordinates from [0,1] to box size for all snapshots or
        for snapshot `timestep`, using the h-matrix to treat orthogonal or
        triclinic boxes.

        Parameters
        ----------
        timestep : int, optional
            The timestep of the snapshot being unscaled.
        """
        if timestep == -1:
            print("Unscaling dump ...")
            x = self.names["x"]
            y = self.names["y"]
            z = self.names["z"]
            for snap in self.snaps:
                self.unscale_one(snap, x, y, z)
        else:
            i = self.find_snapshot(timestep)
            x = self.names["x"]
            y = self.names["y"]
            z = self.names["z"]
            self.unscale_one(self.snaps[i], x, y, z)

    def unscale_one(self, snap: Snapshot, x: int, y: int, z: int) -> None:
        """
        Unscale a single snapshot.

        Parameters
        ----------
        snap : Snapshot
            Snapshot being unscaled.
        x : int
            Column id of the x coordinate.
        y : int
            Column id of the y coordinate.
        z : int
            Column id of the z coordinate.
        """
        if not snap.triclinic:
            xprd = snap.xhi - snap.xlo
            yprd = snap.yhi - snap.ylo
            zprd = snap.zhi - snap.zlo
            atoms = snap.atoms
            if atoms is not None:
                atoms[:, x] = snap.xlo + atoms[:, x]*xprd
                atoms[:, y] = snap.ylo + atoms[:, y]*yprd
                atoms[:, z] = snap.zlo + atoms[:, z]*zprd
        else:
            # Creating h-matrix:
            h_mat, _ = self.unit_cell(snap)
            atoms = snap.atoms
            if atoms is not None:
                atoms[:, x] = snap.xlo + atoms[:, x] * \
                    h_mat[0] + atoms[:, y]*h_mat[5] + atoms[:, z]*h_mat[4]
                atoms[:, y] = snap.ylo + atoms[:, y]*h_mat[1] + \
                    atoms[:, z]*h_mat[3]
                atoms[:, z] = snap.zlo + atoms[:, z]*h_mat[2]

    def wrap(self) -> None:
        """
        Wraps coordinates from outside box to inside.
        """
        print("Wrapping dump ...")
        x = self.names["x"]
        y = self.names["y"]
        z = self.names["z"]
        ix = self.names["ix"]
        iy = self.names["iy"]
        iz = self.names["iz"]
        for snap in self.snaps:
            xprd = snap.xhi - snap.xlo
            yprd = snap.yhi - snap.ylo
            zprd = snap.zhi - snap.zlo
            atoms = snap.atoms
            atoms[:, x] -= atoms[:, ix]*xprd
            atoms[:, y] -= atoms[:, iy]*yprd
            atoms[:, z] -= atoms[:, iz]*zprd

    def unwrap(self) -> None:
        """
        Unwraps coordinates from inside box to outside.
        """
        print("Unwrapping dump ...")
        x = self.names["x"]
        y = self.names["y"]
        z = self.names["z"]
        ix = self.names["ix"]
        iy = self.names["iy"]
        iz = self.names["iz"]
        for snap in self.snaps:
            xprd = snap.xhi - snap.xlo
            yprd = snap.yhi - snap.ylo
            zprd = snap.zhi - snap.zlo
            atoms = snap.atoms
            atoms[:, x] += atoms[:, ix]*xprd
            atoms[:, y] += atoms[:, iy]*yprd
            atoms[:, z] += atoms[:, iz]*zprd

    @staticmethod
    def write_header(snap: Snapshot, file: IO, cols: str) -> None:
        """
        Write the 'header' section of a dump file.

        Parameters
        ----------
        snap : Snapshot
            The snapshot for which header is written
        file : TextIO
            Pointer to the file to which the header is written
        cols : str
            The space-separated column names
        """
        file.write("ITEM: TIMESTEP\n")
        file.write(str(snap.time) + "\n")
        file.write("ITEM: NUMBER OF ATOMS\n")
        file.write(str(snap.n_atoms) + "\n")
        if snap.boxstr:
            file.write(snap.boxstr)
        else:
            file.write("ITEM: BOX BOUNDS\n")
        if snap.triclinic:
            file.write(f"{snap.xlo} {snap.xhi} {snap.xy}\n")
            file.write(f"{snap.ylo} {snap.yhi} {snap.xz}\n")
            file.write(f"{snap.zlo} {snap.zhi} {snap.yz}\n")
        else:
            file.write(f"{snap.xlo} {snap.xhi}\n")
            file.write(f"{snap.ylo} {snap.yhi}\n")
            file.write(f"{snap.zlo} {snap.zhi}\n")
        file.write("ITEM: ATOMS " + cols + "\n")

    def write(self, filename: str, mode: str = "w"
              ) -> None:
        """
        write a single dump file from current selection.

        Parameters
        ----------
        filename: str
            _description_
        mode : str, optional
            Whether write to a new file or append to an existing file.
        """
        col_str = self.names2str()
        if "id" in self.names:
            atom_id_col = self.names["id"]
        else:
            raise ValueError("atom id column not found!")
        if "type" in self.names:
            atom_type_col = self.names["type"]
        else:
            raise ValueError("atom type column not found!")
        with open(filename, mode) as snapshot:
            for snap in self.snaps:
                self.write_header(snap, snapshot, col_str)
                for atom in snap.atoms:
                    line = ""
                    for j in range(snap.n_props):
                        if j == atom_id_col or j == atom_type_col:
                            line += str(int(atom[j])) + " "
                        else:
                            line += str(atom[j]) + " "
                    snapshot.write(line)
                    snapshot.write("\n")
        print(f"{self.n_snaps} snapshots are written to a single dump file.")

    def scatter(self, prefix: str) -> None:
        """
        Write one dump file per snapshot from the current selection, using
        `prefix` as the prefix for filenames.

        Parameters
        ----------
        prefix: str
            The prefix used in the filenames.
        """
        col_str = self.names2str()
        for snap in self.snaps:
            print(snap.time, flush=True)
            filename = prefix + "." + str(snap.time) + '.lammpstrj'
            with open(filename, "w") as snapshot:
                self.write_header(snap, snapshot, col_str)
                for atom in snap.atoms:
                    line = ""
                    for j in range(snap.n_props):
                        if j < 2:
                            line += str(int(atom[j])) + " "
                        else:
                            line += str(atom[j]) + " "
                    snapshot.write(line)
        print(f"{self.n_snaps} snapshot(s) are written to "
              f"{self.n_snaps} file(s).")
