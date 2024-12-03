import glob
import warnings
import os
from typing import Callable, List, Tuple, Optional, Union, Dict, Any

import numpy as np

from ..manage.typer import ParserInstance, GroupT, PropertyT
from ..manage.utilizer import invalid_keyword, sort_filenames
from ..manage import organizer
from .distributions import distributions_generator
from .correlations import acf_generator


class AnalyzerBase:
    """
    Base class for analyzing simulation data with various methods for time-series,
    histograms, and other statistical analyses.

    Attributes
    ----------
    parser : ParserType
        A parser to process filenames and infer metadata.
    group : str
        Group type ('bug', 'nucleoid', 'all').
    geometry : str
        The geometry of the simulation box (e.g., 'cylindrical', 'cubic').
    topology : str
        The topology of the polymer.
    is_segment : bool
        Whether the input data is in segments or as wholes.
    methods : Dict[str, Callable]
        A dictionary mapping property types to their corresponding methods.
    """

    def __init__(
        self,
        input_database: str,
        hierarchy: str,
        has_stamp: bool,
        properties: Dict[str, Tuple[Any, ...]],
        group: GroupT,
        is_segment: bool,
    ):
        self._input_database = input_database
        self._hierarchy = hierarchy
        self._has_stamp = has_stamp
        self._group = group
        self._is_segment = is_segment
        self._properties = properties
        self.methods = self._initialize_methods()

    @property
    def parser(self) -> ParserInstance:
        """
        Returns the parser object containing metadata about the simulation.
        """
        if self._parser is None:
            raise AttributeError("'_parser' has not been initialized.")
        return self._parser

    @property
    def input_database(self) -> str:
        """
        Path to the input_database; a 'space' directory at a given 'phase'.
        """
        return self._input_database

    @property
    def hierarchy(self) -> str:
        """
        Hierarchy of the directories and files within the `input_database`;
        for instance, "/N*/N*" means files that starts with "N" and are
        located in directories starting with "N".
        """
        return self._hierarchy

    @property
    def has_stamp(self) -> bool:
        """
        Whether `artifacts` have 'stamp' files (True) or 'whole' (False).
        """
        return self._has_stamp

    @property
    def group(self) -> GroupT:
        """
        Return the current group.
        """
        return self._group

    @property
    def is_segment(self) -> bool:
        """
        Whether `artifacts` are 'segment' (True) or 'whole' (False).
        """
        return self._is_segment

    def _initialize_methods(self) -> Dict[str, Callable]:
        """
        Initializes the dictionary mapping property types to their corresponding methods.

        Returns
        -------
        methods : dict
            Dictionary where keys are property types, and values are methods.
        """
        return {
            "tseries_properties": self.time_series,
            "acf_tseries_properties": self.time_series,
            "hist_properties": self.histograms,
            "hist_properties_no_edge": self.histograms,
            "hist2d_properties": self.histograms_2d,
            "hist2d_edges": self.histograms_2d,
            "rho_phi_hist_properties": self.histograms,
            "nonscalar_hist_t_properties": self.nonscalar_time_series,
            "nonscalar_mat_t_properties": self.nonscalar_time_series,
        }

    def analyze_measures(self) -> None:
        """read in the 'probe' artifacts of the 'group' particles based on the
        `hierarchy` of directories and files from the `input_database` path to the
        'probe' phase of a 'space' and creates the 'analysis' phase at that parent
        directory of the 'probe' of that 'space', infers 'space' names from
        `input_database` path and creates a 'space' directories at various stages
        in the 'analysis' directory for both 'bug' and 'all' groups.

        `tseries_properties`, `hists_properties`, `rho_phi_hists_properties` are
        list of tuples in which each tuple has three string members. The first
        string is the name of a physical property, the second one is the particle
        type, and the last one is `group` type.

        Parameters
        ----------
        input_database: str
            Path to the input_database; a 'space' directory at a given 'phase'.
        hierarchy: str
            Hierarchy of the directories and files within the `input_database`;
            for instance, "/N*/N*" means files that starts with "N" and are
            located in directories starting with "N".
        parser: ParserT
            A class from 'PolyPhys.manage.parser' module that parses filenames
            or filepaths to infer information about a file.
        group : {'bug', 'nucleoid', 'all'}, default cylindrical
            Shape of the simulation box.
        geometry : {'cylindrical', 'slit', 'cubic'}
            Shape of the simulation box.
        topology:
            Topology of the polymer
        is_segment: bool
            Whether `artifacts` are 'segment' (True) or 'whole' (False)
        has_stamp: bool.
            Whether `artifacts` have 'stamp' files (True) or 'whole' (False).
        """
        artifacts = glob.glob(self.input_database + self.hierarchy)
        if not artifacts:
            raise ValueError(
                f"No files found in '{self.input_database + self.hierarchy}'."
            )
        self.artifacts = artifacts
        # Define save paths
        self.save_to_ens = organizer.make_database(
            self.input_database, 'analysis', stage='ens', group=self.group
        )
        self.save_to_ens_avg = organizer.make_database(
            self.input_database, 'analysis', stage='ensAvg', group=self.group
        )
        self.save_to_whole = None
        if self.is_segment is True: 
            self.save_to_whole = organizer.make_database(
                self.input_database, 'analysis', stage='wholeSim', group=self.group
            )
        # Handle stamp files if present
        if self.has_stamp is True:
            self._stamps()

        # Dynamically call methods based on property types in kwargs
        for prop_type, prop_value in kwargs.items():
            if prop_value is not None:
                method = self.methods.get(prop_type)
                if method:
                    method(**kwargs)

    def time_series(
        self,
        properties: List[Tuple[str, str, str]],
    ) -> None:
        """Runs various statistical analyses on `artifacts` of
        each of `TimeSeriesT` types in a given `geometry` and then
        writes the ensembles and ensemble-averages of time series and its
        associated analyses to file.

        If the `is_segment` is `True`, `artifacts` are "segments" and
        are "vertically" merged to create "wholes".

        Issue
        -----
        If the order of the acf_generator and the first ensemble_avg is
        change, the column names in each ensemble dataframe in 'ensembles'
        changes.

        Parameters
        ----------
        properties: list of HistogramT
            A list of tuples where each tuple has three members: the property
            name, entity, and group of a 'time-series' property.
        """
        for prop_, entity, group in properties:
            tseries = sort_filenames(
                self.artifacts,
                fmts=['-' + prop_ + entity + '.npy']
            )
            # if the type is 'segment' we need to merge files and create
            # 'whole' files.
            if self.is_segment is True:
                wholes = organizer.whole_from_segments(
                    prop_ + entity,
                    tseries,  # type: ignore
                    self.parser,
                    group,
                    relation='tseries',
                    save_to=self.save_to_whole
                )
            else:
                wholes = organizer.whole_from_file(
                    tseries,  # type: ignore
                    self.parser,
                    group,
                )
            ensembles = organizer.ensemble(
                    prop_ + entity,
                    wholes,
                    self.parser,
                    group,
                    'vector',
                    save_to=self.save_to_ens
            )
            _ = organizer.ensemble_avg(
                    prop_ + entity,
                    ensembles,
                    group,
                    'dataframe',
                    save_to=self.save_to_ens_avg
            )

    def time_series_afc(
        self,
        properties: List[Tuple[str, str, str]],
    ) -> None:
        """Runs various statistical analyses on `artifacts` of
        each of `TimeSeriesT` types in a given `geometry` and then
        writes the ensembles and ensemble-averages of time series and its
        associated analyses to file.

        If the `is_segment` is `True`, `artifacts` are "segments" and
        are "vertically" merged to create "wholes".

        Issue
        -----
        If the order of the acf_generator and the first ensemble_avg is
        change, the column names in each ensemble dataframe in 'ensembles'
        changes.

        Parameters
        ----------
        properties: list of HistogramT
            A list of tuples where each tuple has three members: the property
            name, entity, and group of a 'time-series' property. For
            `cls_tseries_properties`, the auto correlation function (AFC) is
            also computed.
        n_lags: int, default 7000
            Maximum lag in the auto correlation function (AFC).
        alpha: float, default 0.05
            If a number is given, the confidence intervals for the given level
            are returned. For instance if alpha=.05, 95 % confidence intervals
            are returned where the standard deviation is computed according to
            Bartlett's formula.
        """
        for prop_, entity, group, n_lags, alpha in properties:
            tseries = sort_filenames(
                self.artifacts,
                fmts=['-' + prop_ + entity + '.npy']
            )
            # if the type is 'segment' we need to merge files and create
            # 'whole' files.
            if self.is_segment is True:
                wholes = organizer.whole_from_segments(
                    prop_ + entity,
                    tseries,  # type: ignore
                    self.parser,
                    group,
                    'tseries',
                    save_to=self.save_to_whole
                )
            else:
                wholes = organizer.whole_from_file(
                    tseries,  # type: ignore
                    self.parser,
                    group,
                )
            ensembles = organizer.ensemble(
                    prop_ + entity,
                    wholes,
                    self.parser,
                    group,
                    'vector',
                    save_to=self.save_to_ens
                )
            acfs, lower_cls, upper_cls = acf_generator(
                prop_ + entity,
                ensembles,  # type: ignore
                n_lags,
                alpha,
                group,
                save_to=self.save_to_ens
            )
            _ = organizer.ensemble_avg(
                    prop_ + entity,
                    ensembles,
                    group,
                    'dataframe',
                    save_to=self.save_to_ens_avg
                )
            _ = organizer.ensemble_avg(
                    prop_ + entity + '-acf',
                    acfs,
                    group,
                    'dataframe',
                    save_to=self.save_to_ens_avg
                )
            _ = organizer.ensemble_avg(
                    prop_ + entity + '-acfLowerCi',
                    lower_cls,
                    group,
                    'dataframe',
                    save_to=self.save_to_ens_avg
                )
            _ = organizer.ensemble_avg(
                    prop_ + entity + '-acfUpperCi',
                    upper_cls,
                    group,
                    'dataframe',
                    save_to=self.save_to_ens_avg
                )

    def histogram(
        self,
        properties: List[Tuple[str, str, str]],
    ) -> None:
        """Runs various statistical analyses on `artifacts` of
        each of `HistogramT` types in a given `geometry` and then
        writes the ensembles and ensemble-averages of time series and its
        associated analyses to file.

        If the `segment` is `True`, `artifacts` are "segments" and
        are "horizontally" merged to create "wholes".

        Issue
        -----
        HThis function only work for spatial distributions created by
        SpatialDistribution class.

        Parameters
        ----------
        artifacts: list of str
            List of path to different artifacts generated by 'probed'
            trajectory files.
        parser: ParserT
            A class from 'PolyPhys.manage.parser' module that parses filenames
            or filepaths to infer information about a file.
        geometry : {'cylindrical', 'slit', 'cubic'}, default cylindrical
            Shape of the simulation box.
        topology:
            Topology of the polymer
        is_segment: bool, default False
            Whether `artifacts` are 'segment' or 'whole'
        save_to: tuple of three str
            Absolute or relative path of the directories to which wholes,
            ensembles, and ensemble-averages are saved.
        rho_phi_hist_properties: list of HistogramT, default None
            A list of tuples where each tuple has four members: the direction,
            direction long name, entity, and group of a 'histogram' property.
            These histogram properties are then used to calculate the local
            number density and volume fraction.
        hist_properties: list of HistogramT default None
            A list of tuples where each tuple has three members: the direction,
            entity, and group of a 'histogram' property.
        hist_properties_no_edge: list of HistogramT default None
            A list of tuples where each tuple has three members: the direction,
            entity, and group of a 'histogram' property. This type of histrograms
            does not have an accopanying edge.
        """
        for direction, entity, group in properties:
            hists = sort_filenames(
                self.artifacts,
                fmts=['-' + direction + 'Hist' + entity + '.npy']
            )
            edges = sort_filenames(
                self.artifacts,
                fmts=['-' + direction + 'Edge' + entity + '.npy']
            )
            if self.is_segment is True:
                wholes = organizer.whole_from_segments(
                    direction + 'Hist' + entity,
                    hists,  # type: ignore
                    self.parser,
                    group,
                    'histogram',
                    save_to=self.save_to_whole
                )
                whole_edges = organizer.whole_from_segments(
                    direction + 'Edge' + entity,
                    edges,  # type: ignore
                    self.parser,
                    group,
                    'bin_edge',
                    save_to=self.save_to_whole
                )
            else:
                wholes = organizer.whole_from_file(
                    hists,  # type: ignore
                    self.parser,
                    group,
                )
                whole_edges = organizer.whole_from_file(
                    edges,  # type: ignore
                    self.parser,
                    group,
                )
            ensembles = organizer.ensemble(
                direction + 'Hist' + entity,
                wholes,
                self.parser,
                group,
                'vector',
                whole_edges=whole_edges,
                save_to=self.save_to_ens
            )
            _ = organizer.ensemble_avg(
                direction + 'Hist' + entity,
                ensembles,
                group,
                'dataframe',
                save_to=self.save_to_ens_avg
            )

    def histogram_rho_phi(
        self,
        properties: List[Tuple[str, str, str]],
    ) -> None:
        """Runs various statistical analyses on `artifacts` of
        each of `HistogramT` types in a given `geometry` and then
        writes the ensembles and ensemble-averages of time series and its
        associated analyses to file.

        If the `segment` is `True`, `artifacts` are "segments" and
        are "horizontally" merged to create "wholes".

        Issue
        -----
        HThis function only work for spatial distributions created by
        SpatialDistribution class.

        Parameters
        ----------
        artifacts: list of str
            List of path to different artifacts generated by 'probed'
            trajectory files.
        parser: ParserT
            A class from 'PolyPhys.manage.parser' module that parses filenames
            or filepaths to infer information about a file.
        geometry : {'cylindrical', 'slit', 'cubic'}, default cylindrical
            Shape of the simulation box.
        topology:
            Topology of the polymer
        is_segment: bool, default False
            Whether `artifacts` are 'segment' or 'whole'
        save_to: tuple of three str
            Absolute or relative path of the directories to which wholes,
            ensembles, and ensemble-averages are saved.
        rho_phi_hist_properties: list of HistogramT, default None
            A list of tuples where each tuple has four members: the direction,
            direction long name, entity, and group of a 'histogram' property.
            These histogram properties are then used to calculate the local
            number density and volume fraction.
        hist_properties: list of HistogramT default None
            A list of tuples where each tuple has three members: the direction,
            entity, and group of a 'histogram' property.
        hist_properties_no_edge: list of HistogramT default None
            A list of tuples where each tuple has three members: the direction,
            entity, and group of a 'histogram' property. This type of histrograms
            does not have an accopanying edge.
        """
        for direction, entity, group in properties:
            hists = sort_filenames(
                self.artifacts,
                fmts=['-' + direction + 'Hist' + entity + '.npy']
            )
            edges = sort_filenames(
                self.artifacts,
                fmts=['-' + direction + 'Edge' + entity + '.npy']
            )
            if self.is_segment is True:
                wholes = organizer.whole_from_segments(
                    direction + 'Hist' + entity,
                    hists,
                    self.parser,
                    group,
                    'histogram',
                    save_to=self.save_to_whole
                )
                whole_edges = organizer.whole_from_segments(
                    direction + 'Edge' + entity,
                    edges,
                    self.parser,
                    group,
                    'bin_edge',
                    save_to=self.save_to_whole
                )
            else:
                wholes = organizer.whole_from_file(
                    hists,  # type: ignore
                    self.parser,
                    group,
                )
                whole_edges = organizer.whole_from_file(
                    edges,  # type: ignore
                    self.parser,
                    group,
                )
            # 'whole' dataframes, each with a 'whole' columns
            rho_wholes, phi_wholes = distributions_generator(
                wholes,
                whole_edges,
                group,
                entity,
                direction,
                self.parser,  # type: ignore
                save_to=self.save_to_whole
            )
            ensembles = organizer.ensemble(
                direction + 'Hist' + entity,
                wholes,
                self.parser,
                group,
                'vector',
                whole_edges=whole_edges,
                save_to=self.save_to_ens
            )
            _ = organizer.ensemble_avg(
                direction + 'Hist' + entity,
                ensembles,
                group,
                'dataframe',
                save_to=self.save_to_ens_avg
            )
            ensembles = organizer.ensemble(
                direction + 'Rho' + entity,
                rho_wholes,
                self.parser,
                group,
                'vector',
                whole_edges=whole_edges,
                save_to=self.save_to_ens
            )
            _ = organizer.ensemble_avg(
                direction + 'Rho' + entity,
                ensembles,
                group,
                'dataframe',
                save_to=self.save_to_ens_avg
            )
            ensembles = organizer.ensemble(
                direction + 'Phi' + entity,
                phi_wholes,
                self.parser,
                group,
                'vector',
                whole_edges=whole_edges,
                save_to=self.save_to_ens
            )
            _ = organizer.ensemble_avg(
                direction + 'Phi' + entity,
                ensembles,
                group,
                'dataframe',
                save_to=self.save_to_ens_avg
            )

    def histogram_no_edge(
        self,
        properties: List[Tuple[str, str, str]],
    ) -> None:
        """Runs various statistical analyses on `artifacts` of
        each of `HistogramT` types in a given `geometry` and then
        writes the ensembles and ensemble-averages of time series and its
        associated analyses to file.

        If the `segment` is `True`, `artifacts` are "segments" and
        are "horizontally" merged to create "wholes".

        Issue
        -----
        HThis function only work for spatial distributions created by
        SpatialDistribution class.

        Parameters
        ----------
        artifacts: list of str
            List of path to different artifacts generated by 'probed'
            trajectory files.
        parser: ParserT
            A class from 'PolyPhys.manage.parser' module that parses filenames
            or filepaths to infer information about a file.
        geometry : {'cylindrical', 'slit', 'cubic'}, default cylindrical
            Shape of the simulation box.
        topology:
            Topology of the polymer
        is_segment: bool, default False
            Whether `artifacts` are 'segment' or 'whole'
        save_to: tuple of three str
            Absolute or relative path of the directories to which wholes,
            ensembles, and ensemble-averages are saved.
        rho_phi_hist_properties: list of HistogramT, default None
            A list of tuples where each tuple has four members: the direction,
            direction long name, entity, and group of a 'histogram' property.
            These histogram properties are then used to calculate the local
            number density and volume fraction.
        hist_properties: list of HistogramT default None
            A list of tuples where each tuple has three members: the direction,
            entity, and group of a 'histogram' property.
        hist_properties_no_edge: list of HistogramT default None
            A list of tuples where each tuple has three members: the direction,
            entity, and group of a 'histogram' property. This type of histrograms
            does not have an accopanying edge.
        """
        for direction, entity, group in properties:
            hists = sort_filenames(
                self.artifacts,
                fmts=['-' + direction + 'Hist' + entity + '.npy']
            )
            if self.is_segment is True:
                wholes = organizer.whole_from_segments(
                    direction + 'Hist' + entity,
                    hists,  # type: ignore
                    self.parser,
                    group,
                    'histogram',
                    save_to=self.save_to_whole
                )
            else:
                wholes = organizer.whole_from_file(
                    hists,  # type: ignore
                    self.parser,
                    group,
                )
            ensembles = organizer.ensemble(
                direction + 'Hist' + entity,
                wholes,
                self.parser,
                group,
                'vector',
                save_to=self.save_to_ens
            )
            _ = organizer.ensemble_avg(
                direction + 'Hist' + entity,
                ensembles,
                group,
                'dataframe',
                save_to=self.save_to_ens_avg
            )

    def histogram_2d(
        self,
        properties: List[Tuple[str, str, str]],
    ) -> None:
        """Runs various statistical analyses on `artifacts` of
        each of `HistogramT` types in a given `geometry` and then
        writes the ensembles and ensemble-averages of time series and its
        associated analyses to file.

        If the `segment` is `True`, `artifacts` are "segments" and
        are "horizontally" merged to create "wholes".

        Issue
        -----
        HThis function only work for spatial distributions created by
        SpatialDistribution class.

        Parameters
        ----------
        artifacts: list of str
            List of path to different artifacts generated by 'probed'
            trajectory files.
        parser: ParserT
            A class from 'PolyPhys.manage.parser' module that parses filenames
            or filepaths to infer information about a file.
        geometry : {'cylindrical', 'slit', 'cubic'}, default cylindrical
            Shape of the simulation box.
        topology:
            Topology of the polymer
        is_segment: bool, default False
            Whether `artifacts` are 'segment' or 'whole'
        save_to : tuple of three str
            Absolute or relative path of the directories to which wholes,
            ensembles, and ensemble-averages are saved.
        hist2d_properties: list of HistogramT default None
            A list of tuples where each tuple has three members: the direction,
            entity, and group of a 'histogram' property.
        hist2d_edges: list of HistogramT default None
            A list of tuples where each tuple has three members: the direction,
            entity, and group of a 'edge' property.
        """
        for direction, entity, group in properties:
            hists = sort_filenames(
                self.artifacts,
                fmts=[direction + 'Hist' + entity + '.npy']
            )
            if self.is_segment is True:
                wholes = organizer.whole_from_segments(
                    direction + 'Hist' + entity,
                    hists,  # type: ignore
                    self.parser,
                    group,
                    'histogram',
                    save_to=self.save_to_whole
                )
            else:
                wholes = organizer.whole_from_file(
                    hists,  # type: ignore
                    self.parser,
                    group
                )
            ensembles = organizer.ensemble(
                direction + 'Hist' + entity,
                wholes,
                self.parser,
                group,
                'matrix',
                save_to=self.save_to_ens
            )
            _ = organizer.ensemble_avg(
                direction + 'Hist' + entity,
                ensembles,
                group,
                'ndarray',
                save_to=self.save_to_ens_avg
            )

    def histogram_2d_edges(
        self,
        properties: List[Tuple[str, str]],
    ) -> None:
        """Runs various statistical analyses on `artifacts` of
        each of `HistogramT` types in a given `geometry` and then
        writes the ensembles and ensemble-averages of time series and its
        associated analyses to file.

        If the `segment` is `True`, `artifacts` are "segments" and
        are "horizontally" merged to create "wholes".

        Issue
        -----
        HThis function only work for spatial distributions created by
        SpatialDistribution class.

        Parameters
        ----------
        artifacts: list of str
            List of path to different artifacts generated by 'probed'
            trajectory files.
        parser: ParserT
            A class from 'PolyPhys.manage.parser' module that parses filenames
            or filepaths to infer information about a file.
        geometry : {'cylindrical', 'slit', 'cubic'}, default cylindrical
            Shape of the simulation box.
        topology:
            Topology of the polymer
        is_segment: bool, default False
            Whether `artifacts` are 'segment' or 'whole'
        save_to : tuple of three str
            Absolute or relative path of the directories to which wholes,
            ensembles, and ensemble-averages are saved.
        hist2d_properties: list of HistogramT default None
            A list of tuples where each tuple has three members: the direction,
            entity, and group of a 'histogram' property.
        hist2d_edges: list of HistogramT default None
            A list of tuples where each tuple has three members: the direction,
            entity, and group of a 'edge' property.
        """
        for direction, group in properties:
            edges = sort_filenames(
                self.artifacts,
                fmts=[direction + 'Edge' + '.npy']
            )
            if self.is_segment is True:
                edge_wholes = organizer.whole_from_segments(
                    direction + 'Edge',
                    edges,  # type: ignore
                    self.parser,
                    group,
                    'bin_edge',
                    save_to=self.save_to_whole
                )
            else:
                edge_wholes = organizer.whole_from_file(
                    edges,  # type: ignore
                    self.parser,
                    group
                )
            ensembles = organizer.ensemble(
                direction + 'Edge',
                edge_wholes,
                self.parser,
                group,
                'bin_edge',
                save_to=self.save_to_ens
            )
            _ = organizer.ensemble_avg(
                direction + 'Edge',
                ensembles,
                group,
                'bin_edge',
                save_to=self.save_to_ens_avg
            )

    def nonscalar_hist_t(self, properties) -> None:
        """Runs overs all 'segment' `artifacts` of each of
        "NonScalarTimeSeriesT" types in a given 'geometry', takes time average over
        a given axis of each "NonScalarTimeSeriesT", and then writes the ensembles
        and ensemble-averages of time-averaged  "NonScalarTimeSeriesT" and its
        associated analyses to file. `nonscalar_hist_properties` are time-varying
        histograms while `nonscalar_matrix_properties` are time-varying matrices.

        A non_scalar property is either a time-varying 1D array (a vector) or a
        time-varying 2D one (a matrix). See the "Notes" and "Issues" below.

        Notes
        -----
        Currently, the vector-like properties are histograms collected over time so
        there are passed by `nonscalar_hist_properties` argument. The "avg_axis"
        ,defined below , sets the axis showing the time evolution, thus allowing
        performing time-averaging over a vector and finding the averages of its
        components.

        If the `is_segment` is `True`, `artifacts` are "segments" and
        are "vertically" merged to create "wholes".

        Issues
        ------
        1. Based on the Notes above, different type of nonscalar properties should
        be passed to this function and new functions should be defined for
        matrix-like nonscalar properties, so it can be possible to define the
        ensemble and ensemble-averaged versions.

        2. `nonscalar_hist_properties` is currently list of time-varying histogram
        that do not need any more process as those passed to "histograms".

        Parameters
        ----------
        artifacts: list of str
            List of path to different artifacts generated by 'probed'
            trajectory files.
        parser: ParserT
            A class from 'PolyPhys.manage.parser' module that parses filenames
            or filepaths to infer information about a file.
        geometry : {'cylindrical', 'slit', 'cubic'}
            Shape of the simulation box.
        topology:
            Topology of the polymer
        is_segment: bool
            Whether `artifacts` are 'segment' or 'whole'
        save_to : tuple of three str
            Absolute or relative path of the directories to which wholes,
            ensembles, and ensemble-averages are saved.
        nonscalar_hist_t_properties: list of NonScalarTimeSeriesT, default None
            A list of tuples in which each tuple has four string members. The
            first string is the name of a physical property, the second one is
            the particle type, the third one is `group` type, and the last one
            is the axis over which these physical properties are all of
            nonscalar form.
        nonscalar_mat_t_properties: list of NonScalarTimeSeriesT, default None
            A list of tuples in which each tuple has three string members. The
            first string is the name of a physical property, the second one is
            the particle type, and the last one is `group` type. These physical
            properties are all of nonscalar form.
        """
        for prop_, entity, group, avg_axis in properties:
            tseries = sort_filenames(
                    self.artifacts,
                    fmts=['-' + prop_ + entity + '.npy']
                )
            if self.is_segment is True:
                wholes = organizer.whole_from_segments(
                    prop_ + entity,
                    tseries,  # type: ignore
                    self.parser,
                    group,
                    'tseries',
                    save_to=self.save_to_whole
                )
            else:
                wholes = organizer.whole_from_file(
                    tseries,  # type: ignore
                    self.parser,
                    group,  
                )
            # Cleaning process on "clustersHistTFoci":
            # See `polyphys.analyze.clusters.count_clusters` for why we have to
            # drop the first item and what it means.
            # See `polyphys.probe.prober.trans_fuci_bug` for how
            # `count_clusters` applied.
            # Time-averaging process:
            wholes = {
                whole_name: np.mean(whole_array, axis=avg_axis)
                for whole_name, whole_array in wholes.items()
            }
            # changing prop_ name after averaging:
            prop_old = prop_
            prop_ = ''.join(prop_.split('T'))  # type: ignore
            warnings.warn(
                f"property name '{prop_old}' changed to"
                f" '{prop_}' after averaging over time.",
                UserWarning
            )
            ensembles = organizer.ensemble(
                    prop_ + entity,
                    wholes,
                    self.parser,
                    group,
                    'vector',
                    save_to=self.save_to_ens
            )
            _ = organizer.ensemble_avg(
                    prop_ + entity,
                    ensembles,
                    group,
                    'dataframe',
                    save_to=self.save_to_ens_avg
            )

    def nonscalar_cluster_hist_t(self, properties) -> None:
        """Runs overs all 'segment' `artifacts` of each of
        "NonScalarTimeSeriesT" types in a given 'geometry', takes time average over
        a given axis of each "NonScalarTimeSeriesT", and then writes the ensembles
        and ensemble-averages of time-averaged  "NonScalarTimeSeriesT" and its
        associated analyses to file. `nonscalar_hist_properties` are time-varying
        histograms while `nonscalar_matrix_properties` are time-varying matrices.

        A non_scalar property is either a time-varying 1D array (a vector) or a
        time-varying 2D one (a matrix). See the "Notes" and "Issues" below.

        Notes
        -----
        Currently, the vector-like properties are histograms collected over time so
        there are passed by `nonscalar_hist_properties` argument. The "avg_axis"
        ,defined below , sets the axis showing the time evolution, thus allowing
        performing time-averaging over a vector and finding the averages of its
        components.

        If the `is_segment` is `True`, `artifacts` are "segments" and
        are "vertically" merged to create "wholes".

        Issues
        ------
        1. Based on the Notes above, different type of nonscalar properties should
        be passed to this function and new functions should be defined for
        matrix-like nonscalar properties, so it can be possible to define the
        ensemble and ensemble-averaged versions.

        2. `nonscalar_hist_properties` is currently list of time-varying histogram
        that do not need any more process as those passed to "histograms".

        Parameters
        ----------
        artifacts: list of str
            List of path to different artifacts generated by 'probed'
            trajectory files.
        parser: ParserT
            A class from 'PolyPhys.manage.parser' module that parses filenames
            or filepaths to infer information about a file.
        geometry : {'cylindrical', 'slit', 'cubic'}
            Shape of the simulation box.
        topology:
            Topology of the polymer
        is_segment: bool
            Whether `artifacts` are 'segment' or 'whole'
        save_to : tuple of three str
            Absolute or relative path of the directories to which wholes,
            ensembles, and ensemble-averages are saved.
        nonscalar_hist_t_properties: list of NonScalarTimeSeriesT, default None
            A list of tuples in which each tuple has four string members. The
            first string is the name of a physical property, the second one is
            the particle type, the third one is `group` type, and the last one
            is the axis over which these physical properties are all of
            nonscalar form.
        nonscalar_mat_t_properties: list of NonScalarTimeSeriesT, default None
            A list of tuples in which each tuple has three string members. The
            first string is the name of a physical property, the second one is
            the particle type, and the last one is `group` type. These physical
            properties are all of nonscalar form.
        """
        #cluster_list = \
        #        ['clustersHistT', 'clustersHistDangleT', 'clustersHistBridgeT',
        #        'clustersHistDirDepT']
        for prop_, entity, group, avg_axis in properties:
            tseries = sort_filenames(
                    self.artifacts,
                    fmts=['-' + prop_ + entity + '.npy']
                )
            if self.is_segment is True:
                wholes = organizer.whole_from_segments(
                    prop_ + entity,
                    tseries,  # type: ignore
                    self.parser,
                    group,
                    'tseries',
                    save_to=self.save_to_whole
                )
            else:
                wholes = organizer.whole_from_file(
                    tseries,  # type: ignore
                    self.parser,
                    group,  
                )
            # Cleaning process on "clustersHistTFoci":
            # See `polyphys.analyze.clusters.count_clusters` for why we have to
            # drop the first item and what it means.
            # See `polyphys.probe.prober.trans_fuci_bug` for how
            # `count_clusters` applied.
            
            wholes = {
                whole_name: whole_array[:, 1:]
                for whole_name, whole_array in wholes.items()
            }
            # Time-averaging process:
            wholes = {
                whole_name: np.mean(whole_array, axis=avg_axis)
                for whole_name, whole_array in wholes.items()
            }
            # changing prop_ name after averaging:
            prop_old = prop_
            prop_ = ''.join(prop_.split('T'))  # type: ignore
            warnings.warn(
                f"property name '{prop_old}' changed to"
                f" '{prop_}' after averaging over time.",
                UserWarning
            )
            ensembles = organizer.ensemble(
                    prop_ + entity,
                    wholes,
                    self.parser,
                    group,
                    'vector',
                    save_to=self.save_to_ens
            )
            _ = organizer.ensemble_avg(
                    prop_ + entity,
                    ensembles,
                    group,
                    'dataframe',
                    save_to=self.save_to_ens_avg
            )

    def clusters_time_series(self, properties) -> None:
        """Runs overs all 'segment' `artifacts` of each of
        "NonScalarTimeSeriesT" types in a given 'geometry', takes time average over
        a given axis of each "NonScalarTimeSeriesT", and then writes the ensembles
        and ensemble-averages of time-averaged  "NonScalarTimeSeriesT" and its
        associated analyses to file. `nonscalar_hist_properties` are time-varying
        histograms while `nonscalar_matrix_properties` are time-varying matrices.

        A non_scalar property is either a time-varying 1D array (a vector) or a
        time-varying 2D one (a matrix). See the "Notes" and "Issues" below.

        Parameters
        ----------
        properties: list of NonScalarTimeSeriesT, default None
            A list of tuples in which each tuple has three string members. The
            first string is the name of a physical property, the second one is
            the particle type, and the last one is `group` type. These physical
            properties are all of nonscalar form.
        """
        for prop_, entity, group in properties:
            whole_paths = sort_filenames(
                self.artifacts,
                fmts=['-' + prop_ + entity + '.npy']
            )
            # changing prop_ name after averaging:
            # Here not only "T" is dropped but also the "property name" is
            # fully changed.
            #if prop_ == 'distMatT' and \
            #    parser.__name__ in ['TransFociCub', 'TransFociCyl',
            #                        'SumRuleCubHeteroRing',
            #                        'SumRuleCubHeteroLinear']:  # type: ignore
            wholes_hists, wholes_rdfs, wholes_tseries = \
                organizer.whole_from_dist_mat_t(
                    whole_paths,
                    self.parser,
                    group
                    )
            # Foci hists:
            # "wholes_hists" are of 'dataframe' whole_type, we directly get
            # its ensemble-averaged properties.
            # changing prop_ name after averaging:
            prop_new = 'pairDistHist'
            warnings.warn(
                f"property name '{prop_}' changed to"
                f" '{prop_new}' after averaging over time.",
                UserWarning
            )
            _ = organizer.ensemble(
                    prop_new + entity + '-ensAvg',
                    wholes_hists,
                    self.parser,
                    group,
                    'dataframe',
                    save_to=self.save_to_ens_avg
            )
            # Foci rdfs:
            # "wholes_rdfs" are of 'dataframe' whole_type, we directly get
            # its ensemble-averaged properties.
            # changing prop_ name after averaging:
            prop_new = 'pairDistRdf'
            warnings.warn(
                f"property name '{prop_}' changed to"
                f" '{prop_new}' after averaging over time.",
                UserWarning
            )
            _ = organizer.ensemble(
                    prop_new + entity + '-ensAvg',
                    wholes_rdfs,
                    self.parser,
                    group,
                    'dataframe',
                    save_to=self.save_to_ens_avg
            )
            # Foci time series:
            # "wholes_tseries" are of 'dataframe' whole_type, we directly
            # get its ensemble-averaged properties.
            # changing prop_ name after averaging:
            prop_new = 'pairDistT'
            warnings.warn(
                f"property name '{prop_}' changed to"
                f" '{prop_new}' after averaging over time.",
                UserWarning
            )
            # For 'dataframe' whole_type, we directly get the
            # ensemble-averaged properties.
            _ = organizer.ensemble(
                    prop_new + entity + '-ensAvg',
                    wholes_tseries,
                    self.parser,
                    group,
                    'dataframe',
                    save_to=self.save_to_ens_avg
            )

    def nonscalar_time_series(self, properties) -> None:
        """Runs overs all 'segment' `artifacts` of each of
        "NonScalarTimeSeriesT" types in a given 'geometry', takes time average over
        a given axis of each "NonScalarTimeSeriesT", and then writes the ensembles
        and ensemble-averages of time-averaged  "NonScalarTimeSeriesT" and its
        associated analyses to file. `nonscalar_hist_properties` are time-varying
        histograms while `nonscalar_matrix_properties` are time-varying matrices.

        A non_scalar property is either a time-varying 1D array (a vector) or a
        time-varying 2D one (a matrix). See the "Notes" and "Issues" below.

        Parameters
        ----------
        properties: list of NonScalarTimeSeriesT, default None
            A list of tuples in which each tuple has three string members. The
            first string is the name of a physical property, the second one is
            the particle type, and the last one is `group` type. These physical
            properties are all of nonscalar form.
        """
        for prop_, entity, group in properties:
            tseries = sort_filenames(
                self.artifacts, fmts=["-" + prop_ + entity + ".npy"]
            )
            if self.is_segment is True:
                wholes = organizer.whole_from_segments(
                    prop_ + entity,
                    tseries,
                    self.parser,
                    self.group,
                    "tseries",
                    save_to=self.save_to_whole,
                )
            else:
                wholes = organizer.whole_from_file(
                    tseries, self.parser, self.group  # type: ignore
                )
            ensembles = organizer.ensemble(
                prop_ + entity,
                wholes,
                self.parser,
                self.group,
                "matrix",
                save_to=self.save_to_ens,
            )
            _ = organizer.ensemble_avg(
                prop_ + entity,
                ensembles,
                group,
                "ndarray",
                save_to=self.save_to_ens_avg,
            )

    def _stamps(self) -> None:
        """
        Combine stamps files.
        """
        stamp_files = sort_filenames(self.artifacts, fmts=["-stamps.csv"])
        if self.is_segment is True:
            segments_stamps = organizer.children_stamps(
                stamp_files,
                "segment",  # lineage of the children stamps
                save_to=self.save_to_whole,  # save segment stamps
            )
            whole_stamps = organizer.parents_stamps(
                segments_stamps,
                "segment",  # lineage of the children stamps
                save_to=self.save_to_ens,  # save whole stamps
            )
        else:
            whole_stamps = organizer.children_stamps(
                stamp_files,
                "whole",  # lineage of the children stamps
                save_to=self.save_to_ens,  # save whole stamps
            )
        _ = organizer.parents_stamps(
            whole_stamps,
            "whole",  # lineage of the children stamps
            save_to=self.save_to_ens_avg,  # save ensemble-averaged stamps
        )
