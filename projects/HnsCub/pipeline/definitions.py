"""
HnsCub project focuses on the behaviour of a linear
heterogeneous in the presence of crowders in free/unconfined space.
"""
from typing import Dict, List

lineage_attributes: Dict[str, Dict[str, str]] = {
        "segment": {  # lcube twice of l
            "nmon": "N",
            "bend_mm": "kbmm",
            "eps_hc": "epshc",
            "nhns": "nh",
            "dcrowd": "ac",
            "ncrowd": "nc",
            "lcube": "l",
            "ensemble_id": "ens",
            "segment_id": "j",
        },
        "whole": {  # lcube twice of l
            "nmon": "N",
            "bend_mm": "kbmm",
            "eps_hc": "epshc",
            "nhns": "nh",
            "dcrowd": "ac",
            "ncrowd": "nc",
            "lcube": "l",
            "ensemble_id": "ens",
        },
        "ensemble_long": {  # lcube twice of l
            "nmon": "N",
            "bend_mm": "kbmm",
            "eps_hc": "epshc",
            "nhns": "nh",
            "dcrowd": "ac",
            "ncrowd": "nc",
            "lcube": "l",
        },
        "ensemble": {
            "nmon": "N",
            "bend_mm": "kbmm",
            "eps_hm": "epshm",
            "nhns": "nh",
            "dcrowd": "ac",
            "ncrowd": "nc",
        },
        "space": {
            "nmon": "N",
            "bend_mm": "kbmm",
            "eps_hc": "epshc",
            "nhns": "nh",
            "dcrowd": "ac",
        },
    }
physical_attributes: Dict[str, List[str]] = {
        "segment": [
            "dmon", "mmon", "dhns", "mhns", "mcrowd", "eps_others",
            "phi_m_bulk", "rho_m_bulk", "phi_c_bulk", "rho_c_bulk",
            "phi_hns_bulk", "rho_hns_bulk", "dt", "ndump", "adump"
            ],
        "whole": [
            "dmon", "mmon", "dhns", "mhns", "mcrowd", "eps_others",
            "phi_m_bulk", "rho_m_bulk", "phi_c_bulk", "rho_c_bulk",
            "phi_hns_bulk", "rho_hns_bulk",  "dt", "ndump", "adump"
            ],
        "ensemble_long": [
            "dmon", "mmon", "dhns", "mhns", "mcrowd", "eps_others",
            "phi_m_bulk", "rho_m_bulk", "phi_c_bulk", "rho_c_bulk",
            "phi_hns_bulk", "rho_hns_bulk", "dt", "ndump", "adump"
            ],
        "ensemble": [
            "dmon", "mmon", "dhns", "mhns", "mcrowd", "eps_others",
            "dt", "ndump", "adump"
            ],
        "space": [
            "dmon", "mmon", "dhns", "mhns", "mcrowd", "eps_others",
            "dt", "ndump", "adump"
            ]
    }

project_definition = {
    {"stage": "analysis",
     "hierarchy": "al*/al*",  # dir/file
    "parser": SumRuleCubHeteroLinear,
    "group": "bug",
    "geometry": "cubic",
    "topology": "linear",
    "is_segment": False,
    "has_stamp": True,
    "nonscalar_hist_t_properties": [
        # property_, species, group, avg_axis
        ("bondsHistT", "Foci", "bug", 0),
        ("clustersHistT", "Foci", "bug", 0)
    ],
    "nonscalar_mat_t_properties": [
        # property_, species, group, avg_axis
        ("distMatT", "Foci", "bug"),
        ("principalT", "Mon", "bug")
    ],
    "acf_tseries_properties": [
        ("gyrT", "Mon", "bug"),
        ("rfloryT", "Mon", "bug"),
        ("shapeT", "Mon", "bug"),
        ("asphericityT", "Mon", "bug")
    ]
}
ANALYSIS_DETAILS_ALL = {
    "SumRuleCubHeteroLinearWhole": {
        "hierarchy": "al*/al*",  # dir/file
        "parser": SumRuleCubHeteroLinear,
        "group": "all",
        "geometry": "cubic",
        "topology": "linear",
        "is_segment": True,
        "has_stamp": False,
        "rho_phi_hist_properties": [
            # direction, species, group
            ("r", "Crd", "all"),
            ("r", "Mon", "all"),
            ("r", "Foci", "all"),
        ],
        "hist_properties": [
            # direction, species, group
            ("r", "Dna", "all"),
            ("r", "Crd", "all"),
            ("r", "Mon", "all"),
            ("r", "Foci", "all")],
        "hist2d_properties": [
            # direction, species, group
            ("xy", "Crd", "all"),
            ("xy", "Mon", "all"),
            ("xy", "Dna", "all"),
            ("xy", "Foci", "all"),
            ("xz", "Crd", "all"),
            ("xz", "Mon", "all"),
            ("xz", "Dna", "all"),
            ("xz", "Foci", "all"),
            ("yz", "Crd", "all"),
            ("yz", "Mon", "all"),
            ("yz", "Dna", "all"),
            ("yz", "Foci", "all"),
        ],
        "hist2d_edges": [
            # direction, group
            ("x", "all"),
            ("y", "all"),
            ("z", "all")
        ]
    },
    