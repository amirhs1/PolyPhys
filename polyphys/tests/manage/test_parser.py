"""Tests for :mod:`polyphys.manage.parser`."""
import os
import re
from collections import OrderedDict
from operator import attrgetter
from pathlib import Path
from typing import Any, NamedTuple

import pandas as pd
import pytest

from polyphys.manage.parser import (
    HnsCub,
    HnsCyl,
    KeywordSpec,
    ParserBase,
    Species,
    SumRuleCubHeteroLinear,
    SumRuleCubHeteroRing,
    SumRuleCubHeteroRingBiDisperse,
    SumRuleCyl,
    TransFociCub,
    TransFociCyl,
    TwoMonDepCub,
)
from polyphys.manage.utils import (
    number_density_cube,
    number_density_cylinder,
    volume_fraction_cube,
    volume_fraction_cylinder,
)


class ParseCase(NamedTuple):
    """A hand-written artifact name and what its parser must give."""
    parser_class: type[ParserBase]
    lineage: Any
    group: Any
    artifact: str
    name: str
    values: dict[str, Any]
    physical: list[str]


class ValidCase(NamedTuple):
    """A valid ``ensemble_long`` artifact name for one parser class."""
    parser_class: type[ParserBase]
    artifact: str
    lineage: Any
    group: Any


@pytest.fixture(params=[
    pytest.param(
        ParseCase(
            TwoMonDepCub, 'segment', 'bug',
            artifact='am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100'
                     'tdump500ens3.j01.bug',
            name='am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100tdump500'
                 'ens3.j01',
            values=dict(dmon=5.0, nmon=2, dcrowd=1.0, ncrowd=100, lcube=20.0,
                        d_sur=2.0, dt=0.01, bdump=50, adump=100, tdump=500,
                        ensemble_id=3, segment_id=1),
            physical=['phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c']),
        id='TwoMonDepCub-segment'),
    pytest.param(
        ParseCase(
            TwoMonDepCub, 'whole', 'bug',
            artifact='am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100'
                     'tdump500ens3.bug',
            name='am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100tdump500'
                 'ens3',
            values=dict(dmon=5.0, nmon=2, dcrowd=1.0, ncrowd=100, lcube=20.0,
                        d_sur=2.0, dt=0.01, bdump=50, adump=100, tdump=500,
                        ensemble_id=3),
            physical=['phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c']),
        id='TwoMonDepCub-whole'),
    pytest.param(
        ParseCase(
            TwoMonDepCub, 'ensemble_long', 'bug',
            artifact='am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100'
                     'tdump500',
            name='am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100tdump500',
            values=dict(dmon=5.0, nmon=2, dcrowd=1.0, ncrowd=100, lcube=20.0,
                        d_sur=2.0, dt=0.01, bdump=50, adump=100, tdump=500),
            physical=['phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c', 'rho_bulk_c']),
        id='TwoMonDepCub-ensemble_long'),
    pytest.param(
        ParseCase(
            TwoMonDepCub, 'ensemble', 'bug',
            artifact='nm2am5.0ac1.0nc100sd2.0',
            name='nm2am5.0ac1.0nc100sd2.0',
            values=dict(dmon=5.0, nmon=2, dcrowd=1.0, ncrowd=100, d_sur=2.0),
            physical=[]),
        id='TwoMonDepCub-ensemble'),
    pytest.param(
        ParseCase(
            TwoMonDepCub, 'space', 'bug',
            artifact='nm2am5.0ac1.0nc100',
            name='nm2am5.0ac1.0nc100',
            values=dict(dmon=5.0, nmon=2, dcrowd=1.0, ncrowd=100),
            physical=[]),
        id='TwoMonDepCub-space'),
    pytest.param(
        ParseCase(
            SumRuleCyl, 'segment', 'bug',
            artifact='N2000epsilon5.0r10.5lz300.0sig1.0nc117000dt0.005'
                     'bdump5000adump10000ens2.j02.bug',
            name='N2000epsilon5.0r10.5lz300.0sig1.0nc117000dt0.005bdump5000'
                 'adump10000ens2.j02',
            values=dict(nmon=2000, epsilon=5.0, dcyl=20, lcyl=600.0,
                        dcrowd=1.0, ncrowd=117000, dt=0.005, bdump=5000,
                        adump=10000, ensemble_id=2, segment_id=2),
            physical=['dmon', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                      'rho_bulk_c']),
        id='SumRuleCyl-segment'),
    pytest.param(
        ParseCase(
            SumRuleCyl, 'whole', 'bug',
            artifact='N2000epsilon5.0r10.5lz300.0sig1.0nc117000dt0.005'
                     'bdump5000adump10000ens2.bug',
            name='N2000epsilon5.0r10.5lz300.0sig1.0nc117000dt0.005bdump5000'
                 'adump10000ens2',
            values=dict(nmon=2000, epsilon=5.0, dcyl=20, lcyl=600.0,
                        dcrowd=1.0, ncrowd=117000, dt=0.005, bdump=5000,
                        adump=10000, ensemble_id=2),
            physical=['dmon', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                      'rho_bulk_c']),
        id='SumRuleCyl-whole'),
    pytest.param(
        ParseCase(
            SumRuleCyl, 'ensemble_long', 'bug',
            artifact='N2000epsilon5.0r10.5lz300.0sig1.0nc117000dt0.005'
                     'bdump5000adump10000',
            name='N2000epsilon5.0r10.5lz300.0sig1.0nc117000dt0.005bdump5000'
                 'adump10000',
            values=dict(nmon=2000, epsilon=5.0, dcyl=20, lcyl=600.0,
                        dcrowd=1.0, ncrowd=117000, dt=0.005, bdump=5000,
                        adump=10000),
            physical=['dmon', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                      'rho_bulk_c']),
        id='SumRuleCyl-ensemble_long'),
    pytest.param(
        ParseCase(
            SumRuleCyl, 'ensemble', 'bug',
            artifact='N2000D20.0ac1.0nc117000',
            name='N2000D20.0ac1.0nc117000',
            values=dict(dmon=1.0, nmon=2000, dcrowd=1.0, ncrowd=117000,
                        dcyl=20),
            physical=['dmon']),
        id='SumRuleCyl-ensemble'),
    pytest.param(
        ParseCase(
            SumRuleCyl, 'space', 'bug',
            artifact='N2000D20.0ac1.0',
            name='N2000D20.0ac1.0',
            values=dict(dmon=1.0, nmon=2000, dcrowd=1.0, dcyl=20),
            physical=['dmon']),
        id='SumRuleCyl-space'),
    pytest.param(
        ParseCase(
            SumRuleCubHeteroRing, 'segment', 'all',
            artifact='al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000'
                     'ens4.j10.ring.all',
            name='al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000ens4.'
                 'j10.ring',
            values=dict(dmon_large=5.0, nmon_large=5, mmon_large=125.0,
                        nmon_small=400, dcrowd=3.0, ncrowd=7735, lcube=90.0,
                        dt=0.005, bdump=2000, adump=5000, ensemble_id=4,
                        segment_id=10),
            physical=['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                      'rho_bulk_m_small', 'phi_bulk_m_large',
                      'rho_bulk_m_large', 'phi_bulk_m', 'rho_bulk_m',
                      'phi_bulk_c', 'rho_bulk_c']),
        id='SumRuleCubHeteroRing-segment'),
    pytest.param(
        ParseCase(
            SumRuleCubHeteroRing, 'whole', 'all',
            artifact='al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000'
                     'ens4.ring.all',
            name='al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000ens4.'
                 'ring',
            values=dict(dmon_large=5.0, nmon_large=5, mmon_large=125.0,
                        nmon_small=400, dcrowd=3.0, ncrowd=7735, lcube=90.0,
                        dt=0.005, bdump=2000, adump=5000, ensemble_id=4),
            physical=['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                      'rho_bulk_m_small', 'phi_bulk_m_large',
                      'rho_bulk_m_large', 'phi_bulk_m', 'rho_bulk_m',
                      'phi_bulk_c', 'rho_bulk_c']),
        id='SumRuleCubHeteroRing-whole'),
    pytest.param(
        ParseCase(
            SumRuleCubHeteroRing, 'ensemble_long', 'all',
            artifact='al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000',
            name='al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000',
            values=dict(dmon_large=5.0, nmon_large=5, mmon_large=125.0,
                        nmon_small=400, dcrowd=3.0, ncrowd=7735, lcube=90.0,
                        dt=0.005, bdump=2000, adump=5000),
            physical=['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                      'rho_bulk_m_small', 'phi_bulk_m_large',
                      'rho_bulk_m_large', 'phi_bulk_m', 'rho_bulk_m',
                      'phi_bulk_c', 'rho_bulk_c']),
        id='SumRuleCubHeteroRing-ensemble_long'),
    pytest.param(
        ParseCase(
            SumRuleCubHeteroRing, 'ensemble', 'all',
            artifact='ns400nl5al5ac3nc7735',
            name='ns400nl5al5ac3nc7735',
            values=dict(dmon_large=5.0, nmon_large=5, nmon_small=400,
                        dcrowd=3.0, ncrowd=7735),
            physical=['dmon_small', 'mmon_small', 'mcrowd']),
        id='SumRuleCubHeteroRing-ensemble'),
    pytest.param(
        ParseCase(
            SumRuleCubHeteroRing, 'space', 'all',
            artifact='ns400nl5al5ac3',
            name='ns400nl5al5ac3',
            values=dict(dmon_large=5.0, nmon_large=5, nmon_small=400,
                        dcrowd=3.0),
            physical=['dmon_small', 'mmon_small', 'mcrowd']),
        id='SumRuleCubHeteroRing-space'),
    pytest.param(
        ParseCase(
            SumRuleCubHeteroLinear, 'segment', 'all',
            artifact='al6nl5ml216ns400ac3nc7735l45dt0.005bdump2000adump5000'
                     'ens1.j11.linear.all',
            name='al6nl5ml216ns400ac3nc7735l45dt0.005bdump2000adump5000ens1.'
                 'j11.linear',
            values=dict(dmon_large=6.0, nmon_large=5, mmon_large=216.0,
                        nmon_small=400, dcrowd=3.0, ncrowd=7735, lcube=90.0,
                        dt=0.005, bdump=2000, adump=5000, ensemble_id=1,
                        segment_id=11),
            physical=['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                      'rho_bulk_m_small', 'phi_bulk_m_large',
                      'rho_bulk_m_large', 'phi_bulk_m', 'rho_bulk_m',
                      'phi_bulk_c', 'rho_bulk_c']),
        id='SumRuleCubHeteroLinear-segment'),
    pytest.param(
        ParseCase(
            SumRuleCubHeteroLinear, 'whole', 'all',
            artifact='al6nl5ml216ns400ac3nc7735l45dt0.005bdump2000adump5000'
                     'ens1.linear.all',
            name='al6nl5ml216ns400ac3nc7735l45dt0.005bdump2000adump5000ens1.'
                 'linear',
            values=dict(dmon_large=6.0, nmon_large=5, mmon_large=216.0,
                        nmon_small=400, dcrowd=3.0, ncrowd=7735, lcube=90.0,
                        dt=0.005, bdump=2000, adump=5000, ensemble_id=1),
            physical=['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                      'rho_bulk_m_small', 'phi_bulk_m_large',
                      'rho_bulk_m_large', 'phi_bulk_m', 'rho_bulk_m',
                      'phi_bulk_c', 'rho_bulk_c']),
        id='SumRuleCubHeteroLinear-whole'),
    pytest.param(
        ParseCase(
            SumRuleCubHeteroLinear, 'ensemble_long', 'all',
            artifact='al6nl5ml216ns400ac3nc7735l45dt0.005bdump2000adump5000',
            name='al6nl5ml216ns400ac3nc7735l45dt0.005bdump2000adump5000',
            values=dict(dmon_large=6.0, nmon_large=5, mmon_large=216.0,
                        nmon_small=400, dcrowd=3.0, ncrowd=7735, lcube=90.0,
                        dt=0.005, bdump=2000, adump=5000),
            physical=['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                      'rho_bulk_m_small', 'phi_bulk_m_large',
                      'rho_bulk_m_large', 'phi_bulk_m', 'rho_bulk_m',
                      'phi_bulk_c', 'rho_bulk_c']),
        id='SumRuleCubHeteroLinear-ensemble_long'),
    pytest.param(
        ParseCase(
            SumRuleCubHeteroLinear, 'ensemble', 'all',
            artifact='ns400nl5al6ac3nc7735',
            name='ns400nl5al6ac3nc7735',
            values=dict(dmon_large=6.0, nmon_large=5, nmon_small=400,
                        dcrowd=3.0, ncrowd=7735),
            physical=['dmon_small', 'mmon_small', 'mcrowd']),
        id='SumRuleCubHeteroLinear-ensemble'),
    pytest.param(
        ParseCase(
            SumRuleCubHeteroLinear, 'space', 'all',
            artifact='ns400nl5al6ac3',
            name='ns400nl5al6ac3',
            values=dict(dmon_large=6.0, nmon_large=5, nmon_small=400,
                        dcrowd=3.0),
            physical=['dmon_small', 'mmon_small', 'mcrowd']),
        id='SumRuleCubHeteroLinear-space'),
    pytest.param(
        ParseCase(
            TransFociCyl, 'segment', 'bug',
            artifact='epss5epsl5r10.5al3nl5ml27ns400ac1nc8340lz69.5dt0.005'
                     'bdump2000adump5000ens8.j05.ring.bug',
            name='epss5epsl5r10.5al3nl5ml27ns400ac1nc8340lz69.5dt0.005'
                 'bdump2000adump5000ens8.j05.ring',
            values=dict(epsilon_small=5.0, epsilon_large=5.0, dcyl=20.0,
                        dmon_large=3.0, nmon_large=5, mmon_large=27.0,
                        nmon_small=400, dcrowd=1.0, ncrowd=8340, lcyl=139.0,
                        dt=0.005, bdump=2000, adump=5000, ensemble_id=8,
                        segment_id=5),
            physical=['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                      'rho_bulk_m_small', 'phi_bulk_m_large',
                      'rho_bulk_m_large', 'phi_bulk_m', 'rho_bulk_m',
                      'phi_bulk_c', 'rho_bulk_c']),
        id='TransFociCyl-segment'),
    pytest.param(
        ParseCase(
            TransFociCyl, 'whole', 'bug',
            artifact='epss5epsl5r10.5al3nl5ml27ns400ac1nc8340lz69.5dt0.005'
                     'bdump2000adump5000ens8.ring.bug',
            name='epss5epsl5r10.5al3nl5ml27ns400ac1nc8340lz69.5dt0.005'
                 'bdump2000adump5000ens8.ring',
            values=dict(epsilon_small=5.0, epsilon_large=5.0, dcyl=20.0,
                        dmon_large=3.0, nmon_large=5, mmon_large=27.0,
                        nmon_small=400, dcrowd=1.0, ncrowd=8340, lcyl=139.0,
                        dt=0.005, bdump=2000, adump=5000, ensemble_id=8),
            physical=['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                      'rho_bulk_m_small', 'phi_bulk_m_large',
                      'rho_bulk_m_large', 'phi_bulk_m', 'rho_bulk_m',
                      'phi_bulk_c', 'rho_bulk_c']),
        id='TransFociCyl-whole'),
    pytest.param(
        ParseCase(
            TransFociCyl, 'ensemble_long', 'bug',
            artifact='epss5epsl5r10.5al3nl5ml27ns400ac1nc8340lz69.5dt0.005'
                     'bdump2000adump5000',
            name='epss5epsl5r10.5al3nl5ml27ns400ac1nc8340lz69.5dt0.005'
                 'bdump2000adump5000',
            values=dict(epsilon_small=5.0, epsilon_large=5.0, dcyl=20.0,
                        dmon_large=3.0, nmon_large=5, mmon_large=27.0,
                        nmon_small=400, dcrowd=1.0, ncrowd=8340, lcyl=139.0,
                        dt=0.005, bdump=2000, adump=5000),
            physical=['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                      'rho_bulk_m_small', 'phi_bulk_m_large',
                      'rho_bulk_m_large', 'phi_bulk_m', 'rho_bulk_m',
                      'phi_bulk_c', 'rho_bulk_c']),
        id='TransFociCyl-ensemble_long'),
    pytest.param(
        ParseCase(
            TransFociCyl, 'ensemble', 'bug',
            artifact='ns400nl5al3D20.0ac1nc8340',
            name='ns400nl5al3D20.0ac1nc8340',
            values=dict(dcyl=20.0, dmon_large=3.0, nmon_large=5,
                        nmon_small=400, dcrowd=1.0, ncrowd=8340),
            physical=['dmon_small', 'mmon_small', 'mcrowd']),
        id='TransFociCyl-ensemble'),
    pytest.param(
        ParseCase(
            TransFociCyl, 'space', 'bug',
            artifact='ns400nl5al3D20.0ac1',
            name='ns400nl5al3D20.0ac1',
            values=dict(dcyl=20.0, dmon_large=3.0, nmon_large=5,
                        nmon_small=400, dcrowd=1.0),
            physical=['dmon_small', 'mmon_small', 'mcrowd']),
        id='TransFociCyl-space'),
    pytest.param(
        ParseCase(
            TransFociCub, 'segment', 'all',
            artifact='al5nl5ml125ns400ac1nc154699l30dt0.005bdump2000adump5000'
                     'ens5.j13.ring.all',
            name='al5nl5ml125ns400ac1nc154699l30dt0.005bdump2000adump5000'
                 'ens5.j13.ring',
            values=dict(dmon_large=5.0, nmon_large=5, mmon_large=125.0,
                        nmon_small=400, dcrowd=1.0, ncrowd=154699, lcube=60.0,
                        dt=0.005, bdump=2000, adump=5000, ensemble_id=5,
                        segment_id=13),
            physical=['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                      'rho_bulk_m_small', 'phi_bulk_m_large',
                      'rho_bulk_m_large', 'phi_bulk_m', 'rho_bulk_m',
                      'phi_bulk_c', 'rho_bulk_c']),
        id='TransFociCub-segment'),
    pytest.param(
        ParseCase(
            TransFociCub, 'whole', 'all',
            artifact='al5nl5ml125ns400ac1nc154699l30dt0.005bdump2000adump5000'
                     'ens5.ring.all',
            name='al5nl5ml125ns400ac1nc154699l30dt0.005bdump2000adump5000'
                 'ens5.ring',
            values=dict(dmon_large=5.0, nmon_large=5, mmon_large=125.0,
                        nmon_small=400, dcrowd=1.0, ncrowd=154699, lcube=60.0,
                        dt=0.005, bdump=2000, adump=5000, ensemble_id=5),
            physical=['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                      'rho_bulk_m_small', 'phi_bulk_m_large',
                      'rho_bulk_m_large', 'phi_bulk_m', 'rho_bulk_m',
                      'phi_bulk_c', 'rho_bulk_c']),
        id='TransFociCub-whole'),
    pytest.param(
        ParseCase(
            TransFociCub, 'ensemble_long', 'all',
            artifact='al5nl5ml125ns400ac1nc154699l30dt0.005bdump2000adump5000',
            name='al5nl5ml125ns400ac1nc154699l30dt0.005bdump2000adump5000',
            values=dict(dmon_large=5.0, nmon_large=5, mmon_large=125.0,
                        nmon_small=400, dcrowd=1.0, ncrowd=154699, lcube=60.0,
                        dt=0.005, bdump=2000, adump=5000),
            physical=['dmon_small', 'mmon_small', 'mcrowd', 'phi_bulk_m_small',
                      'rho_bulk_m_small', 'phi_bulk_m_large',
                      'rho_bulk_m_large', 'phi_bulk_m', 'rho_bulk_m',
                      'phi_bulk_c', 'rho_bulk_c']),
        id='TransFociCub-ensemble_long'),
    pytest.param(
        ParseCase(
            TransFociCub, 'ensemble', 'all',
            artifact='ns400nl5al5ac1nc154699',
            name='ns400nl5al5ac1nc154699',
            values=dict(dmon_large=5.0, nmon_large=5, nmon_small=400,
                        dcrowd=1.0, ncrowd=154699),
            physical=['dmon_small', 'mmon_small', 'mcrowd']),
        id='TransFociCub-ensemble'),
    pytest.param(
        ParseCase(
            TransFociCub, 'space', 'all',
            artifact='ns400nl5al5ac1',
            name='ns400nl5al5ac1',
            values=dict(dmon_large=5.0, nmon_large=5, nmon_small=400,
                        dcrowd=1.0),
            physical=['dmon_small', 'mmon_small', 'mcrowd']),
        id='TransFociCub-space'),
    pytest.param(
        ParseCase(
            HnsCub, 'segment', 'nucleoid',
            artifact='N200kbmm2nh12ac2l25epshc1nc11937ens2.j5.ring.nucleoid',
            name='N200kbmm2nh12ac2l25epshc1nc11937ens2.j5.ring',
            values=dict(nmon=200, bend_mm=2.0, nhns=12, dcrowd=2.0, lcube=50,
                        eps_hc=1.0, ncrowd=11937, ensemble_id=2, segment_id=5),
            physical=['dmon', 'dhns', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                      'rho_bulk_c', 'phi_bulk_hns', 'rho_bulk_hns', 'dt',
                      'ndump', 'adump', 'eps_hm']),
        id='HnsCub-segment'),
    pytest.param(
        ParseCase(
            HnsCub, 'whole', 'nucleoid',
            artifact='N200kbmm2nh12ac2l25epshc1nc11937ens2.ring.nucleoid',
            name='N200kbmm2nh12ac2l25epshc1nc11937ens2.ring',
            values=dict(nmon=200, bend_mm=2.0, nhns=12, dcrowd=2.0, lcube=50,
                        eps_hc=1.0, ncrowd=11937, ensemble_id=2),
            physical=['dmon', 'dhns', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                      'rho_bulk_c', 'phi_bulk_hns', 'rho_bulk_hns', 'dt',
                      'ndump', 'adump', 'eps_hm']),
        id='HnsCub-whole'),
    pytest.param(
        ParseCase(
            HnsCub, 'ensemble_long', 'nucleoid',
            artifact='N200kbmm2nh12ac2l25epshc1nc11937',
            name='N200kbmm2nh12ac2l25epshc1nc11937',
            values=dict(nmon=200, bend_mm=2.0, nhns=12, dcrowd=2.0, lcube=50,
                        eps_hc=1.0, ncrowd=11937),
            physical=['dmon', 'dhns', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                      'rho_bulk_c', 'phi_bulk_hns', 'rho_bulk_hns', 'dt',
                      'ndump', 'adump', 'eps_hm']),
        id='HnsCub-ensemble_long'),
    pytest.param(
        ParseCase(
            HnsCub, 'ensemble', 'nucleoid',
            artifact='N200kbmm2nh12ac2.0epshc1nc11937',
            name='N200kbmm2nh12ac2.0epshc1nc11937',
            values=dict(nmon=200, bend_mm=2.0, nhns=12, dcrowd=2.0, eps_hc=1.0,
                        ncrowd=11937),
            physical=['dmon', 'dhns', 'dt', 'ndump', 'adump', 'eps_hm']),
        id='HnsCub-ensemble'),
    pytest.param(
        ParseCase(
            HnsCub, 'space', 'nucleoid',
            artifact='N200kbmm2nh12ac2.0epshc1',
            name='N200kbmm2nh12ac2.0epshc1',
            values=dict(nmon=200, bend_mm=2.0, nhns=12, dcrowd=2.0,
                        eps_hc=1.0),
            physical=['dmon', 'dhns', 'dt', 'ndump', 'adump', 'eps_hm']),
        id='HnsCub-space'),
    pytest.param(
        ParseCase(
            HnsCyl, 'segment', 'all',
            artifact='N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152ens1.j01.ring.all',
            name='N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152ens1.j01.ring',
            values=dict(nmon=200, bend_mm=2.0, dcyl=8, nhns=4, dcrowd=1.0,
                        lcyl=150, eps_hc=1.0, ncrowd=1152, ensemble_id=1,
                        segment_id=1),
            physical=['dmon', 'dhns', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                      'rho_bulk_c', 'phi_bulk_hns', 'rho_bulk_hns', 'dt',
                      'ndump', 'adump', 'eps_hm']),
        id='HnsCyl-segment'),
    pytest.param(
        ParseCase(
            HnsCyl, 'whole', 'all',
            artifact='N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152ens1.ring.all',
            name='N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152ens1.ring',
            values=dict(nmon=200, bend_mm=2.0, dcyl=8, nhns=4, dcrowd=1.0,
                        lcyl=150, eps_hc=1.0, ncrowd=1152, ensemble_id=1),
            physical=['dmon', 'dhns', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                      'rho_bulk_c', 'phi_bulk_hns', 'rho_bulk_hns', 'dt',
                      'ndump', 'adump', 'eps_hm']),
        id='HnsCyl-whole'),
    pytest.param(
        ParseCase(
            HnsCyl, 'ensemble_long', 'all',
            artifact='N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152',
            name='N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152',
            values=dict(nmon=200, bend_mm=2.0, dcyl=8, nhns=4, dcrowd=1.0,
                        lcyl=150, eps_hc=1.0, ncrowd=1152),
            physical=['dmon', 'dhns', 'phi_bulk_m', 'rho_bulk_m', 'phi_bulk_c',
                      'rho_bulk_c', 'phi_bulk_hns', 'rho_bulk_hns', 'dt',
                      'ndump', 'adump', 'eps_hm']),
        id='HnsCyl-ensemble_long'),
    pytest.param(
        ParseCase(
            HnsCyl, 'ensemble', 'all',
            artifact='N200D8nh4ac1epshc1.0nc1152',
            name='N200D8nh4ac1epshc1.0nc1152',
            values=dict(nmon=200, dcyl=8, nhns=4, dcrowd=1.0, eps_hc=1.0,
                        ncrowd=1152),
            physical=['dmon', 'dhns', 'dt', 'ndump', 'adump', 'eps_hm']),
        id='HnsCyl-ensemble'),
    pytest.param(
        ParseCase(
            HnsCyl, 'space', 'all',
            artifact='N200D8nh4ac1epshc1.0',
            name='N200D8nh4ac1epshc1.0',
            values=dict(nmon=200, dcyl=8, nhns=4, dcrowd=1.0, eps_hc=1.0),
            physical=['dmon', 'dhns', 'dt', 'ndump', 'adump', 'eps_hm']),
        id='HnsCyl-space'),
    pytest.param(
        ParseCase(
            SumRuleCubHeteroRingBiDisperse, 'segment', 'all',
            artifact='al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005'
                     'bdump2000adump5000ens1.j20.ring.all',
            name='al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005bdump2000'
                 'adump5000ens1.j20.ring',
            values=dict(dmon_large=4.0, nmon_large=5, nmon_small=400,
                        dcrowd_one=1.0, dcrowd_two=5.0, pct_one=40.0,
                        ncrowd_one=139229, ncrowd_two=1671, lcube=90.0,
                        dt=0.005, bdump=2000, adump=5000, ensemble_id=1,
                        segment_id=20),
            physical=['dmon_small', 'mmon_small', 'mmon_large', 'mcrowd_one',
                      'mcrowd_two', 'pct_two', 'phi_bulk_m_small',
                      'rho_bulk_m_small', 'phi_bulk_m_large',
                      'rho_bulk_m_large', 'phi_bulk_m', 'rho_bulk_m',
                      'phi_bulk_c_one', 'rho_bulk_c_one', 'phi_bulk_c_two',
                      'rho_bulk_c_two', 'phi_bulk_c', 'rho_bulk_c']),
        id='SumRuleCubHeteroRingBiDisperse-segment'),
    pytest.param(
        ParseCase(
            SumRuleCubHeteroRingBiDisperse, 'whole', 'bug',
            artifact='al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005'
                     'bdump2000adump5000ens1.ring.bug',
            name='al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005bdump2000'
                 'adump5000ens1.ring',
            values=dict(dmon_large=4.0, nmon_large=5, nmon_small=400,
                        dcrowd_one=1.0, dcrowd_two=5.0, pct_one=40.0,
                        ncrowd_one=139229, ncrowd_two=1671, lcube=90.0,
                        dt=0.005, bdump=2000, adump=5000, ensemble_id=1),
            physical=['dmon_small', 'mmon_small', 'mmon_large', 'mcrowd_one',
                      'mcrowd_two', 'pct_two', 'phi_bulk_m_small',
                      'rho_bulk_m_small', 'phi_bulk_m_large',
                      'rho_bulk_m_large', 'phi_bulk_m', 'rho_bulk_m',
                      'phi_bulk_c_one', 'rho_bulk_c_one', 'phi_bulk_c_two',
                      'rho_bulk_c_two', 'phi_bulk_c', 'rho_bulk_c']),
        id='SumRuleCubHeteroRingBiDisperse-whole'),
    pytest.param(
        ParseCase(
            SumRuleCubHeteroRingBiDisperse, 'ensemble_long', 'bug',
            artifact='al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005'
                     'bdump2000adump5000',
            name='al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005bdump2000'
                 'adump5000',
            values=dict(dmon_large=4.0, nmon_large=5, nmon_small=400,
                        dcrowd_one=1.0, dcrowd_two=5.0, pct_one=40.0,
                        ncrowd_one=139229, ncrowd_two=1671, lcube=90.0,
                        dt=0.005, bdump=2000, adump=5000),
            physical=['dmon_small', 'mmon_small', 'mmon_large', 'mcrowd_one',
                      'mcrowd_two', 'pct_two', 'phi_bulk_m_small',
                      'rho_bulk_m_small', 'phi_bulk_m_large',
                      'rho_bulk_m_large', 'phi_bulk_m', 'rho_bulk_m',
                      'phi_bulk_c_one', 'rho_bulk_c_one', 'phi_bulk_c_two',
                      'rho_bulk_c_two', 'phi_bulk_c', 'rho_bulk_c']),
        id='SumRuleCubHeteroRingBiDisperse-ensemble_long'),
    pytest.param(
        ParseCase(
            SumRuleCubHeteroRingBiDisperse, 'ensemble', 'all',
            artifact='al4nl5ns400ac1ab5pctc40nc139229nb1671',
            name='al4nl5ns400ac1ab5pctc40nc139229nb1671',
            values=dict(dmon_large=4.0, nmon_large=5, nmon_small=400,
                        dcrowd_one=1.0, dcrowd_two=5.0, pct_one=40.0,
                        ncrowd_one=139229, ncrowd_two=1671),
            physical=['dmon_small', 'mmon_small', 'pct_two']),
        id='SumRuleCubHeteroRingBiDisperse-ensemble'),
    pytest.param(
        ParseCase(
            SumRuleCubHeteroRingBiDisperse, 'space', 'all',
            artifact='al4nl5ns400ac1ab5pctc40',
            name='al4nl5ns400ac1ab5pctc40',
            values=dict(dmon_large=4.0, nmon_large=5, nmon_small=400,
                        dcrowd_one=1.0, dcrowd_two=5.0, pct_one=40.0),
            physical=['dmon_small', 'mmon_small', 'pct_two']),
        id='SumRuleCubHeteroRingBiDisperse-space'),
])
def parse_case(request: pytest.FixtureRequest) -> ParseCase:
    """
    One hand-written name per parser class and lineage, with the values its
    parser must give: the digits in the name, converted as the class Notes
    state (``lcube = 2 l`` or ``2 hl``, ``lcyl = 2 lz``, ``dcyl = 2 r - 1``).
    """
    case: ParseCase = request.param
    return case


@pytest.fixture
def parsed(parse_case: ParseCase) -> ParserBase:
    """The parser built from the name of `parse_case`."""
    return parse_case.parser_class(parse_case.artifact, parse_case.lineage,
                                   parse_case.group)


@pytest.fixture(params=[
    pytest.param(
        ValidCase(TwoMonDepCub,
                  'am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100'
                  'tdump500', 'ensemble_long', 'bug'),
        id='TwoMonDepCub'),
    pytest.param(
        ValidCase(SumRuleCyl,
                  'N2000epsilon5.0r10.5lz300.0sig1.0nc117000dt0.005'
                  'bdump5000adump10000', 'ensemble_long', 'bug'),
        id='SumRuleCyl'),
    pytest.param(
        ValidCase(SumRuleCubHeteroRing,
                  'al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000',
                  'ensemble_long', 'bug'),
        id='SumRuleCubHeteroRing'),
    pytest.param(
        ValidCase(SumRuleCubHeteroLinear,
                  'al6nl5ml216ns400ac3nc7735l45dt0.005bdump2000adump5000',
                  'ensemble_long', 'bug'),
        id='SumRuleCubHeteroLinear'),
    pytest.param(
        ValidCase(TransFociCyl,
                  'epss5epsl5r10.5al3nl5ml27ns400ac1nc8340lz69.5dt0.005'
                  'bdump2000adump5000', 'ensemble_long', 'bug'),
        id='TransFociCyl'),
    pytest.param(
        ValidCase(TransFociCub,
                  'al5nl5ml125ns400ac1nc154699l30dt0.005bdump2000'
                  'adump5000', 'ensemble_long', 'bug'),
        id='TransFociCub'),
    pytest.param(
        ValidCase(HnsCub, 'N200kbmm2nh12ac2l25epshc1nc11937',
                  'ensemble_long', 'nucleoid'),
        id='HnsCub'),
    pytest.param(
        ValidCase(HnsCyl, 'N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152',
                  'ensemble_long', 'nucleoid'),
        id='HnsCyl'),
    pytest.param(
        ValidCase(SumRuleCubHeteroRingBiDisperse,
                  'al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005'
                  'bdump2000adump5000', 'ensemble_long', 'bug'),
        id='SumRuleCubHeteroRingBiDisperse'),
])
def valid_case(request: pytest.FixtureRequest) -> ValidCase:
    """A valid ``ensemble_long`` name for each parser class."""
    case: ValidCase = request.param
    return case


@pytest.fixture(params=[
    pytest.param((parser_class, lineage),
                 id=f"{parser_class.__name__}-{lineage}")
    for parser_class in (TwoMonDepCub, SumRuleCyl, SumRuleCubHeteroRing,
                         SumRuleCubHeteroLinear, TransFociCyl, TransFociCub,
                         HnsCub, HnsCyl, SumRuleCubHeteroRingBiDisperse)
    for lineage in ('segment', 'whole', 'ensemble_long', 'ensemble', 'space')
])
def declared_case(
    request: pytest.FixtureRequest
) -> tuple[type[ParserBase], Any, str, dict[str, Any]]:
    """
    A name built from a class's declared pattern for one lineage, with 10.5
    as the value of every float keyword and 2 as the value of every integer
    keyword, and the values the parser must store for them. With these
    values every volume fraction stays below 1 in both geometries.
    """
    parser_class, lineage = request.param
    specs = [(keyword, parser_class._keywords[keyword])
             for keyword in parser_class._patterns[lineage]]
    name = ''.join(f"{keyword}{'10.5' if spec.kind is float else '2'}"
                   for keyword, spec in specs)
    expected = {spec.attr: (spec.scale * 10.5 + spec.offset
                            if spec.kind is float else 2)
                for _, spec in specs}
    return parser_class, lineage, name, expected


@pytest.fixture
def chain_declaration() -> dict[str, Any]:
    """The class attributes of a valid parser for a chain in a cube."""
    return {
        '_geometry': 'cubic',
        '_topology': 'linear',
        '_groups': ['bug', 'all'],
        '_dump_keywords': {'bug': 'bdump', 'all': 'adump'},
        '_keywords': {
            'N': KeywordSpec('nmon', int),
            'ac': KeywordSpec('dcrowd', float),
            'nc': KeywordSpec('ncrowd', int),
            'l': KeywordSpec('lcube', float, scale=2.0),
            'ens': KeywordSpec('ensemble_id', int),
            'j': KeywordSpec('segment_id', int),
        },
        '_patterns': {
            'segment': ('N', 'ac', 'nc', 'l', 'ens', 'j'),
            'whole': ('N', 'ac', 'nc', 'l', 'ens'),
            'ensemble_long': ('N', 'ac', 'nc', 'l'),
            'ensemble': ('N', 'ac', 'nc'),
            'space': ('N', 'ac'),
        },
        '_defaults': {'dmon': 1.0, 'bdump': 1000, 'adump': 5000},
        '_species': (Species('m', 'nmon', 'dmon'),
                     Species('c', 'ncrowd', 'dcrowd')),
        '_project_attributes': {
            'segment': ['dmon'],
            'whole': ['dmon'],
            'ensemble_long': ['dmon'],
            'ensemble': ['dmon'],
            'space': ['dmon'],
        },
    }


# Parsing hand-written names


def test_parser_sets_values_from_name(parse_case, parsed):
    """Each lineage's name parses to its digits, converted as documented."""
    actual = attrgetter(*parse_case.values)(parsed)
    # Values are parsed from decimal text; equal to float64 round-off.
    assert actual == pytest.approx(tuple(parse_case.values.values()),
                                   rel=1e-12, abs=1e-12)


def test_parser_reduces_artifact_to_lineage_name(parse_case, parsed):
    """The name drops the group and topology suffixes after the lineage."""
    assert parsed.name == parse_case.name


def test_parser_lists_physical_attributes(parse_case, parsed):
    """Each lineage lists the fixed and computed attributes it sets."""
    assert parsed.physical_attributes == parse_case.physical


def test_parser_names_its_project(parse_case, parsed):
    """The project name is the parser class's name."""
    assert parsed.project_name == parse_case.parser_class.__name__


def test_parser_without_path_reports_no_path(parse_case, parsed):
    """A parser built from a name keeps it and has no filepath or ext."""
    assert (parsed.filename, parsed.filepath, parsed.ext, parsed.ispath,
            parsed.lineage, parsed.group) == (
        parse_case.artifact, None, None, False, parse_case.lineage,
        parse_case.group)


def test_from_name_reports_no_path():
    """A parser built with from_name has no filepath or ext."""
    parser = SumRuleCyl.from_name('N200D10.0ac2.0', 'space', 'bug')

    assert (parser.filepath, parser.ext, parser.ispath) == (None, None, False)


def test_parser_reads_name_from_file_path(parse_case, make_temp_file):
    """A parser built from a file path parses the file's stem."""
    full_path = make_temp_file(parse_case.artifact + '.txt')

    parser = parse_case.parser_class(full_path, parse_case.lineage,
                                     parse_case.group, ispath=True)

    assert (parser.filename, parser.filepath, parser.ext, parser.ispath,
            parser.name) == (
        parse_case.artifact + '.txt', Path(os.path.dirname(full_path)),
        '.txt', True, parse_case.name)


def test_parser_reports_empty_ext_for_file_without_suffix(make_temp_file):
    """A file path with no suffix gives an empty ext, not None."""
    full_path = make_temp_file('N200D10ac2')

    parser = SumRuleCyl(full_path, 'space', 'bug', ispath=True)

    assert (parser.ext, parser.filepath, parser.name) == (
        '', Path(os.path.dirname(full_path)), 'N200D10ac2')


@pytest.mark.parametrize('parser_class, lineage, group, artifact, expected', [
    pytest.param(SumRuleCyl, 'ensemble', 'bug', 'N2000D20.0ac1.0nc117000',
                 ['nmon', 'dcyl', 'dcrowd', 'ncrowd'],
                 id='SumRuleCyl-ensemble'),
    pytest.param(HnsCub, 'space', 'nucleoid', 'N200kbmm2nh12ac2.0epshc1',
                 ['nmon', 'bend_mm', 'nhns', 'dcrowd', 'eps_hc'],
                 id='HnsCub-space'),
])
def test_parser_lists_lineage_attributes_in_pattern_order(
        parser_class, lineage, group, artifact, expected):
    """The parsed attributes are listed in the order the name gives them."""
    parser = parser_class(artifact, lineage, group)

    assert parser.lineage_attributes == expected


def test_name_through_unknown_keyword_raises_value_error():
    """Cutting a name at a keyword it lacks raises ValueError naming it."""
    parser = SumRuleCyl('N2000D20.0ac1.0', 'space', 'bug')

    with pytest.raises(ValueError,
                       match="'zz' keyword not found in 'N2000D20.0ac1.0'"):
        parser._name_through('zz')


def test_parser_lists_each_attribute_once(parsed):
    """No attribute appears twice in `attributes`."""
    assert len(parsed.attributes) == len(set(parsed.attributes))


def test_parser_sets_every_listed_attribute(parsed):
    """Every attribute in `attributes` is set on the instance."""
    assert set(parsed.attributes) <= vars(parsed).keys()


@pytest.mark.parametrize('lineage, artifact, expected', [
    pytest.param('segment',
                 'N2000epsilon5.0r10.5lz300.0sig1.0nc117000dt0.005'
                 'bdump5000adump10000ens2.j02.bug',
                 ['segment', 'whole', 'ensemble_long', 'ensemble', 'space'],
                 id='segment'),
    pytest.param('whole',
                 'N2000epsilon5.0r10.5lz300.0sig1.0nc117000dt0.005'
                 'bdump5000adump10000ens2.bug',
                 ['whole', 'ensemble_long', 'ensemble', 'space'],
                 id='whole'),
    pytest.param('ensemble_long',
                 'N2000epsilon5.0r10.5lz300.0sig1.0nc117000dt0.005'
                 'bdump5000adump10000',
                 ['ensemble_long', 'ensemble', 'space'], id='ensemble_long'),
    pytest.param('ensemble', 'N2000D20.0ac1.0nc117000',
                 ['ensemble', 'space'], id='ensemble'),
    pytest.param('space', 'N2000D20.0ac1.0', ['space'], id='space'),
])
def test_lineage_genealogy_lists_lineage_and_parents(lineage, artifact,
                                                     expected):
    """A lineage's genealogy is itself followed by its parent lineages."""
    parser = SumRuleCyl(artifact, lineage, 'bug')

    assert parser.lineage_genealogy == expected


@pytest.mark.parametrize(
    'parser_class, artifact, group, geometry, topology, groups', [
        pytest.param(TwoMonDepCub, 'nm2am5.0ac1.0nc100', 'bug', 'cubic',
                     'atomic', ['bug', 'all'], id='TwoMonDepCub'),
        pytest.param(SumRuleCyl, 'N2000D20.0ac1.0', 'bug', 'cylindrical',
                     'linear', ['bug', 'all'], id='SumRuleCyl'),
        pytest.param(SumRuleCubHeteroRing, 'ns400nl5al5ac3', 'all', 'cubic',
                     'ring', ['bug', 'all'], id='SumRuleCubHeteroRing'),
        pytest.param(SumRuleCubHeteroLinear, 'ns400nl5al6ac3', 'all',
                     'cubic', 'linear', ['bug', 'all'],
                     id='SumRuleCubHeteroLinear'),
        pytest.param(TransFociCyl, 'ns400nl5al3D20.0ac1', 'bug',
                     'cylindrical', 'ring', ['bug', 'all'],
                     id='TransFociCyl'),
        pytest.param(TransFociCub, 'ns400nl5al5ac1', 'all', 'cubic', 'ring',
                     ['bug', 'all'], id='TransFociCub'),
        pytest.param(HnsCub, 'N200kbmm2nh12ac2.0epshc1', 'nucleoid', 'cubic',
                     'ring', ['nucleoid', 'all'], id='HnsCub'),
        pytest.param(HnsCyl, 'N200D8nh4ac1epshc1.0', 'all', 'cylindrical',
                     'ring', ['nucleoid', 'all'], id='HnsCyl'),
        pytest.param(SumRuleCubHeteroRingBiDisperse,
                     'al4nl5ns400ac1ab5pctc40', 'bug', 'cubic', 'ring',
                     ['bug', 'all'], id='SumRuleCubHeteroRingBiDisperse'),
    ])
def test_parser_reports_geometry_topology_and_groups(
        parser_class, artifact, group, geometry, topology, groups):
    """Each parser reports its project's geometry, topology, and groups."""
    parser = parser_class(artifact, 'space', group)

    assert (parser.geometry, parser.topology, parser.groups) == (
        geometry, topology, groups)


@pytest.mark.parametrize('parser_class, artifact, group, dump_keywords', [
    pytest.param(TwoMonDepCub, 'nm2am5.0ac1.0nc100', 'bug',
                 {'bug': 'bdump', 'all': 'adump'}, id='TwoMonDepCub'),
    pytest.param(SumRuleCyl, 'N2000D20.0ac1.0', 'bug',
                 {'bug': 'bdump', 'all': 'adump'}, id='SumRuleCyl'),
    pytest.param(SumRuleCubHeteroRing, 'ns400nl5al5ac3', 'all',
                 {'bug': 'bdump', 'all': 'adump'},
                 id='SumRuleCubHeteroRing'),
    pytest.param(SumRuleCubHeteroLinear, 'ns400nl5al6ac3', 'all',
                 {'bug': 'bdump', 'all': 'adump'},
                 id='SumRuleCubHeteroLinear'),
    pytest.param(TransFociCyl, 'ns400nl5al3D20.0ac1', 'bug',
                 {'bug': 'bdump', 'all': 'adump'}, id='TransFociCyl'),
    pytest.param(TransFociCub, 'ns400nl5al5ac1', 'all',
                 {'bug': 'bdump', 'all': 'adump'}, id='TransFociCub'),
    pytest.param(HnsCub, 'N200kbmm2nh12ac2.0epshc1', 'nucleoid',
                 {'nucleoid': 'ndump', 'all': 'adump'}, id='HnsCub'),
    pytest.param(HnsCyl, 'N200D8nh4ac1epshc1.0', 'all',
                 {'nucleoid': 'ndump', 'all': 'adump'}, id='HnsCyl'),
    pytest.param(SumRuleCubHeteroRingBiDisperse, 'al4nl5ns400ac1ab5pctc40',
                 'bug', {'bug': 'bdump', 'all': 'adump'},
                 id='SumRuleCubHeteroRingBiDisperse'),
])
def test_parser_maps_each_group_to_its_dump_keyword(
        parser_class, artifact, group, dump_keywords):
    """Each group's dump keyword names the period of its own trajectory."""
    parser = parser_class(artifact, 'space', group)

    assert parser.dump_keywords == dump_keywords


@pytest.mark.parametrize('parser_class, artifact, group, keyword, period', [
    pytest.param(SumRuleCyl,
                 'N200epsilon5.0r10.5lz25.0sig2.0nc1000dt0.005bdump2000'
                 'adump5000ens1', 'bug', 'bdump', 2000, id='SumRuleCyl-bug'),
    pytest.param(SumRuleCyl,
                 'N200epsilon5.0r10.5lz25.0sig2.0nc1000dt0.005bdump2000'
                 'adump5000ens1', 'all', 'adump', 5000, id='SumRuleCyl-all'),
    # HnsCub names carry no dump periods; they are fixed in _defaults.
    pytest.param(HnsCub, 'N200kbmm2.0nh8ac1.0l25.0epshc1.0nc0ens1.ring',
                 'nucleoid', 'ndump', 5000, id='HnsCub-nucleoid'),
    pytest.param(HnsCub, 'N200kbmm2.0nh8ac1.0l25.0epshc1.0nc0ens1.ring',
                 'all', 'adump', 10000, id='HnsCub-all'),
])
def test_dump_keyword_reads_the_group_dump_period(
        parser_class, artifact, group, keyword, period):
    """dump_keyword names the attribute holding the group's dump period."""
    parser = parser_class(artifact, 'whole', group)

    assert (parser.dump_keyword, getattr(parser, parser.dump_keyword)) == (
        keyword, period)


@pytest.mark.parametrize('parser_class, artifact, group, attrs, expected', [
    pytest.param(TwoMonDepCub,
                 'am5.0nm2ac1.nc100hl10.0sd2.0dt0.01bdump50adump100tdump500',
                 'bug', ('dcrowd', 'ncrowd'), (1.0, 100), id='TwoMonDepCub'),
    pytest.param(SumRuleCyl,
                 'N2000epsilon5.r10.5lz300.0sig1.0nc117000dt0.005'
                 'bdump5000adump10000',
                 'bug', ('epsilon', 'dcyl'), (5.0, 20.0), id='SumRuleCyl'),
    pytest.param(SumRuleCubHeteroRing,
                 'al5nl5ml125.ns400ac3nc7735l45dt0.005bdump2000adump5000',
                 'bug', ('mmon_large', 'nmon_small'), (125.0, 400),
                 id='SumRuleCubHeteroRing'),
    pytest.param(SumRuleCubHeteroLinear,
                 'al6nl5ml216.ns400ac3nc7735l45dt0.005bdump2000adump5000',
                 'bug', ('mmon_large', 'nmon_small'), (216.0, 400),
                 id='SumRuleCubHeteroLinear'),
    pytest.param(TransFociCyl,
                 'epss5epsl5r10.5al3.nl5ml27.ns400ac1nc8340lz69.5dt0.005'
                 'bdump2000adump5000',
                 'bug', ('dmon_large', 'mmon_large', 'nmon_small'),
                 (3.0, 27.0, 400), id='TransFociCyl'),
    pytest.param(TransFociCub,
                 'al5nl5ml125.ns400ac1nc154699l30dt0.005bdump2000adump5000',
                 'bug', ('mmon_large', 'nmon_small'), (125.0, 400),
                 id='TransFociCub'),
    pytest.param(HnsCub, 'N200kbmm2nh12ac2l25epshc1.nc11937', 'nucleoid',
                 ('eps_hc', 'ncrowd'), (1.0, 11937), id='HnsCub'),
    pytest.param(HnsCyl, 'N200kbmm2r4.5nh4ac1lz75epshc1.nc1152', 'nucleoid',
                 ('eps_hc', 'ncrowd'), (1.0, 1152), id='HnsCyl'),
    pytest.param(SumRuleCubHeteroRingBiDisperse,
                 'al4nl5ns400ac1ab5.pctc40.nc139229nb1671l45dt0.005'
                 'bdump2000adump5000',
                 'bug', ('dcrowd_two', 'pct_one'), (5.0, 40.0),
                 id='SumRuleCubHeteroRingBiDisperse'),
])
def test_parser_reads_value_with_trailing_dot(parser_class, artifact, group,
                                              attrs, expected):
    """A value written with a trailing dot, such as ``1.``, parses."""
    parser = parser_class(artifact, 'ensemble_long', group)

    # float('1.') is exact.
    assert attrgetter(*attrs)(parser) == pytest.approx(expected, rel=1e-12)


def test_declared_pattern_round_trips(declared_case):
    """
    A name built from a declared pattern parses back to the declared types
    and conversions, so float keywords keep non-integer values (property).
    """
    parser_class, lineage, name, expected = declared_case

    actual = attrgetter(*expected)(
        parser_class(name, lineage, parser_class._groups[0]))

    # One multiplication and one addition in float64.
    assert actual == pytest.approx(tuple(expected.values()), rel=1e-12)
    assert tuple(map(type, actual)) == tuple(map(type, expected.values()))


# Error paths of parsing


@pytest.mark.parametrize('parser_class, artifact, group, keyword, attr', [
    pytest.param(TwoMonDepCub,
                 'am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100xx500',
                 'bug', 'tdump', 'tdump', id='TwoMonDepCub-replaced'),
    pytest.param(SumRuleCyl,
                 'N2000epsilon5.0r10.5lz300.0sig1.0nc117000XX0.005'
                 'bdump5000adump10000',
                 'bug', 'dt', 'dt', id='SumRuleCyl-replaced'),
    pytest.param(SumRuleCubHeteroRing,
                 'al5nl5ml125ns400ac3XX7735l45dt0.005bdump2000adump5000',
                 'bug', 'nc', 'ncrowd', id='SumRuleCubHeteroRing-replaced'),
    pytest.param(SumRuleCubHeteroLinear,
                 'al6nl5ml216XX400ac3nc7735l45dt0.005bdump2000adump5000',
                 'bug', 'ns', 'nmon_small',
                 id='SumRuleCubHeteroLinear-replaced'),
    pytest.param(TransFociCyl,
                 'epss5epsl5r10.5al3nl5xx27ns400ac1nc8340lz69.5dt0.005'
                 'bdump2000adump5000',
                 'bug', 'ml', 'mmon_large', id='TransFociCyl-replaced'),
    pytest.param(TransFociCub,
                 'al5nl5ml125ns400ac1nc154699XX30dt0.005bdump2000adump5000',
                 'bug', 'l', 'lcube', id='TransFociCub-replaced'),
    pytest.param(HnsCub, 'N200kbmm2XX12ac2l25epshc1nc11937', 'nucleoid',
                 'nh', 'nhns', id='HnsCub-replaced'),
    pytest.param(HnsCyl, 'N200kbmm2r4.5nh4ac1lz75XX1.0nc1152', 'nucleoid',
                 'epshc', 'eps_hc', id='HnsCyl-replaced'),
    pytest.param(TwoMonDepCub,
                 '5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100tdump500',
                 'bug', 'am', 'dmon', id='TwoMonDepCub-heading-digit'),
    pytest.param(SumRuleCyl,
                 '2000epsilon5.0r10.5lz300.0sig1.0nc117000dt0.005'
                 'bdump5000adump10000',
                 'bug', 'N', 'nmon', id='SumRuleCyl-heading-digit'),
    pytest.param(SumRuleCubHeteroRing,
                 '5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000',
                 'bug', 'al', 'dmon_large',
                 id='SumRuleCubHeteroRing-heading-digit'),
    pytest.param(SumRuleCubHeteroLinear,
                 '6nl5ml216ns400ac3nc7735l45dt0.005bdump2000adump5000',
                 'bug', 'al', 'dmon_large',
                 id='SumRuleCubHeteroLinear-heading-digit'),
    pytest.param(TransFociCyl,
                 '5epsl5r10.5al3nl5ml27ns400ac1nc8340lz69.5dt0.005'
                 'bdump2000adump5000',
                 'bug', 'epss', 'epsilon_small',
                 id='TransFociCyl-heading-digit'),
    pytest.param(TransFociCub,
                 '5nl5ml125ns400ac1nc154699l30dt0.005bdump2000adump5000',
                 'bug', 'al', 'dmon_large', id='TransFociCub-heading-digit'),
    pytest.param(HnsCub, '200kbmm2nh12ac2l25epshc1nc11937', 'nucleoid',
                 'N', 'nmon', id='HnsCub-heading-digit'),
    pytest.param(HnsCyl, '200kbmm2r4.5nh4ac1lz75epshc1.0nc1152', 'nucleoid',
                 'N', 'nmon', id='HnsCyl-heading-digit'),
    pytest.param(SumRuleCubHeteroRingBiDisperse,
                 'al4nl5ns400ac1ab5pctc40nc139229XX1671l45dt0.005'
                 'bdump2000adump5000',
                 'bug', 'nb', 'ncrowd_two',
                 id='SumRuleCubHeteroRingBiDisperse-replaced'),
])
def test_missing_keyword_raises_value_error(parser_class, artifact, group,
                                            keyword, attr):
    """A name without a pattern keyword raises ValueError naming it."""
    message = (f"keyword '{keyword}' for attribute '{attr}' not found in "
               f"'{artifact}' (lineage 'ensemble_long')")

    with pytest.raises(ValueError, match=re.escape(message)):
        parser_class(artifact, 'ensemble_long', group)


@pytest.mark.parametrize('parser_class, artifact, group, keyword', [
    pytest.param(TwoMonDepCub,
                 'am.nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100tdump500',
                 'bug', 'am', id='TwoMonDepCub'),
    pytest.param(SumRuleCyl,
                 'N.epsilon5.0r10.5lz300.0sig1.0nc117000dt0.005'
                 'bdump5000adump10000',
                 'bug', 'N', id='SumRuleCyl'),
    pytest.param(SumRuleCubHeteroRing,
                 'al.nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000',
                 'bug', 'al', id='SumRuleCubHeteroRing'),
    pytest.param(SumRuleCubHeteroLinear,
                 'al.nl5ml216ns400ac3nc7735l45dt0.005bdump2000adump5000',
                 'bug', 'al', id='SumRuleCubHeteroLinear'),
    pytest.param(TransFociCyl,
                 'epss.epsl5r10.5al3nl5ml27ns400ac1nc8340lz69.5dt0.005'
                 'bdump2000adump5000',
                 'bug', 'epss', id='TransFociCyl'),
    pytest.param(TransFociCub,
                 'al5nl.ml125ns400ac1nc154699l30dt0.005bdump2000adump5000',
                 'bug', 'nl', id='TransFociCub'),
    pytest.param(HnsCub, 'N.kbmm2nh12ac2l25epshc1nc11937', 'nucleoid', 'N',
                 id='HnsCub'),
    pytest.param(HnsCyl, 'N.kbmm2r4.5nh4ac1lz75epshc1.0nc1152', 'nucleoid',
                 'N', id='HnsCyl'),
    pytest.param(SumRuleCubHeteroRingBiDisperse,
                 'al4nl5ns400ac1ab5pctc.nc139229nb1671l45dt0.005'
                 'bdump2000adump5000',
                 'bug', 'pctc', id='SumRuleCubHeteroRingBiDisperse'),
])
def test_non_numeric_value_raises_value_error(parser_class, artifact, group,
                                              keyword):
    """A keyword whose value is not a number raises ValueError naming it."""
    message = (f"value of keyword '{keyword}' in '{artifact}' "
               "(lineage 'ensemble_long'): '.' is not a number")

    with pytest.raises(ValueError, match=re.escape(message)):
        parser_class(artifact, 'ensemble_long', group)


@pytest.mark.parametrize('parser_class, artifact, group, keyword, token', [
    pytest.param(HnsCub, 'N200kbmm2nh12ac2l25epshc1nc1.5', 'nucleoid', 'nc',
                 '1.5', id='HnsCub-nc'),
    pytest.param(SumRuleCyl,
                 'N2.5epsilon5.0r10.5lz300.0sig1.0nc117000dt0.005'
                 'bdump5000adump10000',
                 'bug', 'N', '2.5', id='SumRuleCyl-N'),
])
def test_non_integral_integer_raises_value_error(parser_class, artifact,
                                                 group, keyword, token):
    """A non-integral value of an integer keyword raises ValueError."""
    message = (f"value of keyword '{keyword}' in '{artifact}' "
               f"(lineage 'ensemble_long'): '{token}' is not an integer")

    with pytest.raises(ValueError, match=re.escape(message)):
        parser_class(artifact, 'ensemble_long', group)


def test_empty_artifact_raises_value_error(valid_case):
    """An empty artifact name raises ValueError."""
    with pytest.raises(ValueError,
                       match="'artifact' cannot be an empty string."):
        valid_case.parser_class('', valid_case.lineage, valid_case.group)


def test_invalid_group_raises_value_error(valid_case):
    """A group the project does not have raises ValueError."""
    with pytest.raises(ValueError, match='is an invalid option.'):
        valid_case.parser_class(valid_case.artifact, valid_case.lineage,
                                'invalid')


def test_invalid_lineage_raises_value_error(valid_case):
    """A lineage outside the five lineages raises ValueError."""
    with pytest.raises(ValueError, match='is an invalid option.'):
        valid_case.parser_class(valid_case.artifact, 'invalid',
                                valid_case.group)


# Text representations


def test_str_summarizes_artifact(valid_case):
    """`str` lists the name, geometry, group, lineage, topology, project."""
    parser = valid_case.parser_class(valid_case.artifact, valid_case.lineage,
                                     valid_case.group)

    assert str(parser) == (
        f"Artifact:\n"
        f"    Name: '{valid_case.artifact}',\n"
        f"    Geometry: '{parser.geometry}',\n"
        f"    Group: '{valid_case.group}',\n"
        f"    Lineage: '{valid_case.lineage}',\n"
        f"    Topology: '{parser.topology}',\n"
        f"    Project: '{valid_case.parser_class.__name__}'")


def test_repr_names_artifact(valid_case):
    """`repr` names the file, geometry, group, lineage, topology, project."""
    parser = valid_case.parser_class(valid_case.artifact, valid_case.lineage,
                                     valid_case.group)

    assert repr(parser) == (
        f"Artifact(filename='{valid_case.artifact}', "
        f"geometry='{parser.geometry}', group='{valid_case.group}', "
        f"lineage='{valid_case.lineage}', topology='{parser.topology}', "
        f"project='{valid_case.parser_class.__name__}')")


# Bulk values


@pytest.mark.parametrize('parser_class, artifact, suffix, number, size', [
    pytest.param(TwoMonDepCub,
                 'am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100tdump500',
                 'm', 'nmon', 'dmon', id='TwoMonDepCub-m'),
    pytest.param(TwoMonDepCub,
                 'am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100tdump500',
                 'c', 'ncrowd', 'dcrowd', id='TwoMonDepCub-c'),
    pytest.param(SumRuleCubHeteroRing,
                 'al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000',
                 'm_small', 'nmon_small', 'dmon_small',
                 id='SumRuleCubHeteroRing-m_small'),
    pytest.param(SumRuleCubHeteroRing,
                 'al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000',
                 'm_large', 'nmon_large', 'dmon_large',
                 id='SumRuleCubHeteroRing-m_large'),
    pytest.param(SumRuleCubHeteroRing,
                 'al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000',
                 'c', 'ncrowd', 'dcrowd', id='SumRuleCubHeteroRing-c'),
    pytest.param(SumRuleCubHeteroLinear,
                 'al6nl5ml216ns400ac3nc7735l45dt0.005bdump2000adump5000',
                 'm_small', 'nmon_small', 'dmon_small',
                 id='SumRuleCubHeteroLinear-m_small'),
    pytest.param(SumRuleCubHeteroLinear,
                 'al6nl5ml216ns400ac3nc7735l45dt0.005bdump2000adump5000',
                 'm_large', 'nmon_large', 'dmon_large',
                 id='SumRuleCubHeteroLinear-m_large'),
    pytest.param(SumRuleCubHeteroLinear,
                 'al6nl5ml216ns400ac3nc7735l45dt0.005bdump2000adump5000',
                 'c', 'ncrowd', 'dcrowd', id='SumRuleCubHeteroLinear-c'),
    pytest.param(TransFociCub,
                 'al5nl5ml125ns400ac1nc154699l30dt0.005bdump2000adump5000',
                 'm_small', 'nmon_small', 'dmon_small',
                 id='TransFociCub-m_small'),
    pytest.param(TransFociCub,
                 'al5nl5ml125ns400ac1nc154699l30dt0.005bdump2000adump5000',
                 'm_large', 'nmon_large', 'dmon_large',
                 id='TransFociCub-m_large'),
    pytest.param(TransFociCub,
                 'al5nl5ml125ns400ac1nc154699l30dt0.005bdump2000adump5000',
                 'c', 'ncrowd', 'dcrowd', id='TransFociCub-c'),
    pytest.param(HnsCub, 'N200kbmm2nh12ac2l25epshc1nc11937', 'm', 'nmon',
                 'dmon', id='HnsCub-m'),
    pytest.param(HnsCub, 'N200kbmm2nh12ac2l25epshc1nc11937', 'hns', 'nhns',
                 'dhns', id='HnsCub-hns'),
    pytest.param(HnsCub, 'N200kbmm2nh12ac2l25epshc1nc11937', 'c', 'ncrowd',
                 'dcrowd', id='HnsCub-c'),
    pytest.param(SumRuleCubHeteroRingBiDisperse,
                 'al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005'
                 'bdump2000adump5000',
                 'm_small', 'nmon_small', 'dmon_small',
                 id='SumRuleCubHeteroRingBiDisperse-m_small'),
    pytest.param(SumRuleCubHeteroRingBiDisperse,
                 'al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005'
                 'bdump2000adump5000',
                 'm_large', 'nmon_large', 'dmon_large',
                 id='SumRuleCubHeteroRingBiDisperse-m_large'),
    pytest.param(SumRuleCubHeteroRingBiDisperse,
                 'al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005'
                 'bdump2000adump5000',
                 'c_one', 'ncrowd_one', 'dcrowd_one',
                 id='SumRuleCubHeteroRingBiDisperse-c_one'),
    pytest.param(SumRuleCubHeteroRingBiDisperse,
                 'al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005'
                 'bdump2000adump5000',
                 'c_two', 'ncrowd_two', 'dcrowd_two',
                 id='SumRuleCubHeteroRingBiDisperse-c_two'),
])
def test_cube_bulk_values_match_utils(parser_class, artifact, suffix, number,
                                      size):
    """Cubic bulk values are the cube functions of the parsed values."""
    parser = parser_class(artifact, 'ensemble_long', parser_class._groups[0])
    n, d, lcube = attrgetter(number, size, 'lcube')(parser)

    actual = attrgetter(f'rho_bulk_{suffix}', f'phi_bulk_{suffix}')(parser)

    # The same functions on the same inputs; float64 round-off only.
    assert actual == pytest.approx(
        (number_density_cube(n, d, lcube, pbc=True),
         volume_fraction_cube(n, d, lcube, pbc=True)), rel=1e-12)


@pytest.mark.parametrize('parser_class, artifact, suffix, number, size', [
    pytest.param(SumRuleCyl,
                 'N2000epsilon5.0r10.5lz300.0sig1.0nc117000dt0.005'
                 'bdump5000adump10000',
                 'm', 'nmon', 'dmon', id='SumRuleCyl-m'),
    pytest.param(SumRuleCyl,
                 'N2000epsilon5.0r10.5lz300.0sig1.0nc117000dt0.005'
                 'bdump5000adump10000',
                 'c', 'ncrowd', 'dcrowd', id='SumRuleCyl-c'),
    pytest.param(TransFociCyl,
                 'epss5epsl5r10.5al3nl5ml27ns400ac1nc8340lz69.5dt0.005'
                 'bdump2000adump5000',
                 'm_small', 'nmon_small', 'dmon_small',
                 id='TransFociCyl-m_small'),
    pytest.param(TransFociCyl,
                 'epss5epsl5r10.5al3nl5ml27ns400ac1nc8340lz69.5dt0.005'
                 'bdump2000adump5000',
                 'm_large', 'nmon_large', 'dmon_large',
                 id='TransFociCyl-m_large'),
    pytest.param(TransFociCyl,
                 'epss5epsl5r10.5al3nl5ml27ns400ac1nc8340lz69.5dt0.005'
                 'bdump2000adump5000',
                 'c', 'ncrowd', 'dcrowd', id='TransFociCyl-c'),
    pytest.param(HnsCyl, 'N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152', 'm',
                 'nmon', 'dmon', id='HnsCyl-m'),
    pytest.param(HnsCyl, 'N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152', 'hns',
                 'nhns', 'dhns', id='HnsCyl-hns'),
    pytest.param(HnsCyl, 'N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152', 'c',
                 'ncrowd', 'dcrowd', id='HnsCyl-c'),
])
def test_cylinder_bulk_values_match_utils(parser_class, artifact, suffix,
                                          number, size):
    """Cylindrical bulk values are the cylinder functions of the values."""
    parser = parser_class(artifact, 'ensemble_long', parser_class._groups[0])
    n, d, lcyl, dcyl = attrgetter(number, size, 'lcyl', 'dcyl')(parser)

    actual = attrgetter(f'rho_bulk_{suffix}', f'phi_bulk_{suffix}')(parser)

    # The same functions on the same inputs; float64 round-off only.
    assert actual == pytest.approx(
        (number_density_cylinder(n, d, lcyl, dcyl, pbc=True),
         volume_fraction_cylinder(n, d, lcyl, dcyl, pbc=True)), rel=1e-12)


@pytest.mark.parametrize('parser_class, artifact', [
    pytest.param(SumRuleCubHeteroRing,
                 'al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000',
                 id='SumRuleCubHeteroRing'),
    pytest.param(SumRuleCubHeteroLinear,
                 'al6nl5ml216ns400ac3nc7735l45dt0.005bdump2000adump5000',
                 id='SumRuleCubHeteroLinear'),
    pytest.param(TransFociCyl,
                 'epss5epsl5r10.5al3nl5ml27ns400ac1nc8340lz69.5dt0.005'
                 'bdump2000adump5000',
                 id='TransFociCyl'),
    pytest.param(TransFociCub,
                 'al5nl5ml125ns400ac1nc154699l30dt0.005bdump2000adump5000',
                 id='TransFociCub'),
    pytest.param(SumRuleCubHeteroRingBiDisperse,
                 'al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005'
                 'bdump2000adump5000',
                 id='SumRuleCubHeteroRingBiDisperse'),
])
def test_monomer_bulk_values_are_sums_over_sizes(parser_class, artifact):
    """All-monomer bulk values are the sums of the small and large ones."""
    parser = parser_class(artifact, 'ensemble_long', 'bug')
    small = attrgetter('rho_bulk_m_small', 'phi_bulk_m_small')(parser)
    large = attrgetter('rho_bulk_m_large', 'phi_bulk_m_large')(parser)

    actual = attrgetter('rho_bulk_m', 'phi_bulk_m')(parser)

    # One float64 addition each.
    assert actual == pytest.approx(
        (small[0] + large[0], small[1] + large[1]), rel=1e-12)


@pytest.mark.parametrize('parser_class, artifact', [
    pytest.param(SumRuleCubHeteroRing,
                 'al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000',
                 id='SumRuleCubHeteroRing'),
    pytest.param(SumRuleCubHeteroLinear,
                 'al6nl5ml216ns400ac3nc7735l45dt0.005bdump2000adump5000',
                 id='SumRuleCubHeteroLinear'),
    pytest.param(TransFociCyl,
                 'epss5epsl5r10.5al3nl5ml27ns400ac2nc1000lz69.5dt0.005'
                 'bdump2000adump5000',
                 id='TransFociCyl'),
    pytest.param(TransFociCub,
                 'al5nl5ml125ns400ac2nc1000l30dt0.005bdump2000adump5000',
                 id='TransFociCub'),
])
def test_crowder_mass_is_size_cubed(parser_class, artifact):
    """At a long lineage the crowder mass is its size cubed, m = d^3."""
    parser = parser_class(artifact, 'ensemble_long', 'bug')

    # One float64 power.
    assert parser.mcrowd == pytest.approx(parser.dcrowd ** 3, rel=1e-12)


@pytest.mark.parametrize('parser_class, artifact', [
    pytest.param(SumRuleCubHeteroRing, 'ns400nl5al5ac3nc7735',
                 id='SumRuleCubHeteroRing'),
    pytest.param(SumRuleCubHeteroLinear, 'ns400nl5al6ac3nc7735',
                 id='SumRuleCubHeteroLinear'),
    pytest.param(TransFociCyl, 'ns400nl5al3D20.0ac1nc8340',
                 id='TransFociCyl'),
    pytest.param(TransFociCub, 'ns400nl5al5ac1nc154699', id='TransFociCub'),
])
def test_crowder_mass_is_placeholder_at_ensemble(parser_class, artifact):
    """At the ensemble lineage the crowder mass keeps the placeholder -1."""
    parser = parser_class(artifact, 'ensemble', 'bug')

    assert parser.mcrowd == -1


# HnsCub float keywords


def test_hnscub_keeps_fractional_kbmm_and_epshc():
    """Fractional kbmm and epshc stay floats in every lineage name."""
    parser = HnsCub('N200kbmm2.5nh8ac1.0l25.0epshc1.5nc0ens1.ring', 'whole',
                    'nucleoid')

    assert (parser.bend_mm, parser.eps_hc, parser.whole,
            parser.ensemble_long, parser.ensemble, parser.space) == (
        2.5, 1.5, 'N200kbmm2.5nh8ac1.0l25.0epshc1.5nc0ens1',
        'N200kbmm2.5nh8ac1.0l25.0epshc1.5nc0',
        'N200kbmm2.5nh8ac1epshc1.5nc0', 'N200kbmm2.5nh8ac1epshc1.5')


@pytest.mark.parametrize('lineage', [
    pytest.param('whole', id='whole'),
    pytest.param('ensemble_long', id='ensemble_long'),
    pytest.param('ensemble', id='ensemble'),
    pytest.param('space', id='space'),
])
@pytest.mark.parametrize('epshc, other_epshc', [
    pytest.param('1.5', '1.0', id='half-step'),
    # Six significant digits, as the '%g' format prints, would merge these.
    pytest.param('1.000001', '1', id='seventh-digit'),
])
def test_hnscub_names_differing_in_epshc_stay_distinct(lineage, epshc,
                                                       other_epshc):
    """Names that differ only in epshc keep distinct lineage names."""
    parser = HnsCub(f'N200kbmm2.5nh8ac1.0l25.0epshc{epshc}nc0ens1.ring',
                    'whole', 'nucleoid')
    other = HnsCub(f'N200kbmm2.5nh8ac1.0l25.0epshc{other_epshc}nc0ens1.ring',
                   'whole', 'nucleoid')

    assert getattr(parser, lineage) != getattr(other, lineage)


# Lineage-name rules


@pytest.mark.parametrize(
    'parser_class, group, lineage, artifact, ensemble, space',
    [
        pytest.param(HnsCub, 'nucleoid', 'whole',
                     'N200kbmm2nh12ac2l25epshc1nc0ens1.ring',
                     'N200kbmm2nh12ac2epshc1nc0', 'N200kbmm2nh12ac2epshc1',
                     id='HnsCub-integer-text'),
        pytest.param(HnsCub, 'nucleoid', 'whole',
                     'N200kbmm2.0nh12ac2.0l25.0epshc1.0nc0ens1.ring',
                     'N200kbmm2nh12ac2epshc1nc0', 'N200kbmm2nh12ac2epshc1',
                     id='HnsCub-float-text'),
        pytest.param(HnsCub, 'all', 'whole',
                     'N200kbmm2nh12ac2l25epshc0nc7162ens1.ring',
                     'N200kbmm2nh12ac2epshc0nc7162', 'N200kbmm2nh12ac2epshc0',
                     id='HnsCub-zero'),
        pytest.param(HnsCub, 'all', 'whole',
                     'N200kbmm2nh12ac2l25epshc0.0nc7162ens1.ring',
                     'N200kbmm2nh12ac2epshc0nc7162', 'N200kbmm2nh12ac2epshc0',
                     id='HnsCub-zero-float-text'),
        pytest.param(SumRuleCyl, 'all', 'segment',
                     'N200epsilon5.0r10.5lz25.0sig2.0nc1000dt0.005bdump2000'
                     'adump5000ens1.j03',
                     'N200D20ac2nc1000', 'N200D20ac2',
                     id='SumRuleCyl-converted-dcyl'),
    ],
)
def test_short_name_prints_whole_number_float_as_integer(
        parser_class, group, lineage, artifact, ensemble, space):
    """A whole-number float, zero included, prints as an integer."""
    parser = parser_class(artifact, lineage, group)

    assert (parser.ensemble, parser.space) == (ensemble, space)


@pytest.mark.parametrize(
    'parser_class, group, lineage, artifact, ensemble, space',
    [
        pytest.param(HnsCub, 'nucleoid', 'whole',
                     'N200kbmm2.5nh8ac1.5l25epshc0.25nc0ens1.ring',
                     'N200kbmm2.5nh8ac1.5epshc0.25nc0',
                     'N200kbmm2.5nh8ac1.5epshc0.25', id='HnsCub'),
        pytest.param(SumRuleCyl, 'all', 'segment',
                     'N200epsilon5.0r10.75lz25.0sig1.000001nc1000dt0.005'
                     'bdump2000adump5000ens1.j03',
                     'N200D20.5ac1.000001nc1000', 'N200D20.5ac1.000001',
                     id='SumRuleCyl-converted-dcyl'),
    ],
)
def test_short_name_prints_fractional_float_as_repr(
        parser_class, group, lineage, artifact, ensemble, space):
    """A fractional float prints in full, as repr, in short names."""
    parser = parser_class(artifact, lineage, group)

    assert (parser.ensemble, parser.space) == (ensemble, space)


@pytest.mark.parametrize('short_lineage', [
    pytest.param('ensemble', id='ensemble'),
    pytest.param('space', id='space'),
])
def test_short_name_reparses_to_itself(valid_case, short_lineage):
    """A short name parsed at its own lineage rebuilds the same name."""
    parser = valid_case.parser_class(valid_case.artifact, valid_case.lineage,
                                     valid_case.group)
    name = getattr(parser, short_lineage)

    reparsed = valid_case.parser_class(name, short_lineage, valid_case.group)

    assert getattr(reparsed, short_lineage) == name


@pytest.mark.parametrize('parser_class, group, segment, whole, expected', [
    pytest.param(HnsCub, 'all',
                 'N200kbmm2nh12ac2l25epshc0nc7162ens1.j07.ring.all',
                 'N200kbmm2nh12ac2l25epshc0nc7162ens1.ring.all',
                 'N200kbmm2nh12ac2l25epshc0nc7162ens1', id='HnsCub-ring'),
    pytest.param(SumRuleCubHeteroLinear, 'all',
                 'al6nl5ml216ns400ac3nc7735l45dt0.005bdump2000adump5000'
                 'ens1.j02.linear.all',
                 'al6nl5ml216ns400ac3nc7735l45dt0.005bdump2000adump5000'
                 'ens1.linear.all',
                 'al6nl5ml216ns400ac3nc7735l45dt0.005bdump2000adump5000ens1',
                 id='SumRuleCubHeteroLinear-linear'),
])
def test_whole_name_drops_topology_suffix_from_segment_and_whole(
        parser_class, group, segment, whole, expected):
    """A simulation's segment and whole files give one suffix-free whole."""
    from_segment = parser_class(segment, 'segment', group)
    from_whole = parser_class(whole, 'whole', group)

    assert (from_segment.whole, from_whole.whole) == (expected, expected)


def _stamp_mismatches(path: Path) -> list[tuple[str, ...]]:
    """
    Compare the names in one SciOps ``HnsCub`` stamp file with the rebuilt
    names, under the mapping that the Usage page, "Lineage names", documents.

    A stamp file is ``<lineage_name>-<group>-stamps.csv``; its columns start
    with ``lineage_name`` and the lineage names, with a ``segment`` column
    only for a segment. The expected rebuilt name is the stamp's ``segment``
    or ``whole`` without ``.ring``, its ``ensemble_long`` unchanged, and its
    ``ensemble`` or ``space`` parsed at that lineage and rebuilt; a short
    key must also parse to the artifact's own values. Returns one
    ``(file, column, stamp, expected, rebuilt)`` tuple per mismatch.
    """
    group = path.name.removesuffix('-stamps.csv').rsplit('-', 1)[1]
    frame = pd.read_csv(path, dtype=str)
    lineage = 'segment' if 'segment' in frame.columns else 'whole'
    mismatches = []
    for row in frame.to_dict('records'):
        parser = HnsCub(row['lineage_name'], lineage, group)
        expected = {column: row[column].removesuffix('.ring')
                    for column in ('segment', 'whole') if column in row}
        expected['ensemble_long'] = row['ensemble_long']
        for column in ('ensemble', 'space'):
            old = HnsCub(row[column], column, group)
            expected[column] = getattr(old, column)
            values = attrgetter(*HnsCub._genealogy_attributes[column])
            if values(old) != values(parser):
                mismatches.append((path.name, f'{column} values',
                                   row[column], str(values(parser)),
                                   str(values(old))))
        mismatches.extend(
            (path.name, column, row[column], name, getattr(parser, column))
            for column, name in expected.items()
            if getattr(parser, column) != name)
    return mismatches


@pytest.fixture
def write_stamp(tmp_path: Path):
    """
    Return a factory that writes a one-row stamp file in the SciOps format,
    ``<lineage_name>-<group>-stamps.csv``, from a mapping of column to value,
    and returns its path.
    """
    def _write(group: str, columns: dict[str, str]) -> Path:
        path = tmp_path / f"{columns['lineage_name']}-{group}-stamps.csv"
        path.write_text(','.join(columns) + '\n'
                        + ','.join(columns.values()) + '\n')
        return path
    return _write


@pytest.mark.parametrize('group, columns, flagged', [
    pytest.param(
        'nucleoid',
        {'lineage_name': 'N200kbmm2nh12ac2l25epshc0nc7162ens1.ring',
         'whole': 'N200kbmm2nh12ac2l25epshc0nc7162ens1.ring',
         'ensemble_long': 'N200kbmm2nh12ac2l25epshc0nc7162',
         'ensemble': 'N200kbmm2nh12ac2.0epshc0nc7162',
         'space': 'N200kbmm2nh12ac2.0epshc0', 'nmon': '200', 'n_frames': '10'},
        [], id='whole-matches'),
    pytest.param(
        'all',
        {'lineage_name': 'N200kbmm2nh12ac2l25epshc0nc7162ens1.j07.ring',
         'segment': 'N200kbmm2nh12ac2l25epshc0nc7162ens1.j07.ring',
         'whole': 'N200kbmm2nh12ac2l25epshc0nc7162ens1',
         'ensemble_long': 'N200kbmm2nh12ac2l25epshc0nc7162',
         'ensemble': 'N200kbmm2nh12ac2.0epshc0nc7162',
         'space': 'N200kbmm2nh12ac2.0epshc0', 'nmon': '200', 'n_frames': '10'},
        [], id='segment-matches'),
    pytest.param(
        'nucleoid',
        {'lineage_name': 'N200kbmm2nh12ac2l25epshc0nc7162ens1.ring',
         'whole': 'N200kbmm2nh12ac2l25epshc0nc7162ens1.ring',
         'ensemble_long': 'N200kbmm2nh12ac2l25epshc0nc7162',
         'ensemble': 'N200kbmm2nh12ac2.0epshc0nc7162',
         'space': 'N200kbmm2nh12ac2.0epshc1', 'nmon': '200', 'n_frames': '10'},
        ['space values', 'space'], id='space-value-differs'),
])
def test_stamp_mismatches_flags_only_departures_from_mapping(
        write_stamp, group, columns, flagged):
    """A stamp differing from the rebuilt names only by the mapping passes."""
    path = write_stamp(group, columns)

    mismatches = _stamp_mismatches(path)

    assert [mismatch[1] for mismatch in mismatches] == flagged


@pytest.mark.private_data
def test_hnscub_names_match_existing_stamps_under_mapping(  # pragma: no cover
        private_data_dir):
    """Rebuilt HnsCub names equal the SciOps stamp names under the mapping."""
    stamps = sorted(
        path for path in private_data_dir.rglob('*-stamps.csv')
        if any(part.startswith('HnsCub')
               for part in path.relative_to(private_data_dir).parts[:-1]))

    mismatches = [mismatch for path in stamps
                  for mismatch in _stamp_mismatches(path)]

    assert stamps, 'no stamp file under a directory named HnsCub*'
    assert mismatches == []


# SumRuleCubHeteroRingBiDisperse


@pytest.mark.parametrize('artifact, lineage, group, expected', [
    pytest.param('runs/al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005'
                 'bdump2000adump5000ens1.j20.ring.all.lammpstrj',
                 'segment', 'all',
                 'al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005bdump2000'
                 'adump5000ens1.j20.ring', id='segment-trajectory'),
    pytest.param('runs/al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005'
                 'bdump2000adump5000ens1.ring.bug.data',
                 'whole', 'bug',
                 'al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005bdump2000'
                 'adump5000ens1.ring', id='whole-data'),
    pytest.param('runs/al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005'
                 'bdump2000adump5000ens1.ring.bug.restart.lammpstrj',
                 'whole', 'bug',
                 'al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005bdump2000'
                 'adump5000ens1.ring', id='whole-restart'),
])
def test_bidisperse_reduces_file_path_to_lineage_name(artifact, lineage,
                                                      group, expected):
    """A project file path reduces to its name before the group keyword."""
    parser = SumRuleCubHeteroRingBiDisperse(artifact, lineage, group,
                                            ispath=True)

    assert parser.name == expected


def test_bidisperse_segment_names_its_parents():
    """A segment gives every parent name, short ones without a trailing .0."""
    parser = SumRuleCubHeteroRingBiDisperse(
        'al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005bdump2000adump5000'
        'ens1.j20.ring.all', 'segment', 'all')

    assert (parser.segment, parser.whole, parser.ensemble_long,
            parser.ensemble, parser.space) == (
        'al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005bdump2000adump5000'
        'ens1.j20',
        'al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005bdump2000adump5000'
        'ens1',
        'al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005bdump2000adump5000',
        'al4nl5ns400ac1ab5pctc40nc139229nb1671',
        'al4nl5ns400ac1ab5pctc40')


def test_bidisperse_crowder_bulk_values_are_sums_over_species():
    """All-crowder bulk values are the sums of the two species' values."""
    parser = SumRuleCubHeteroRingBiDisperse(
        'al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005bdump2000adump5000',
        'ensemble_long', 'bug')
    one = attrgetter('rho_bulk_c_one', 'phi_bulk_c_one')(parser)
    two = attrgetter('rho_bulk_c_two', 'phi_bulk_c_two')(parser)

    actual = attrgetter('rho_bulk_c', 'phi_bulk_c')(parser)

    # One float64 addition each.
    assert actual == pytest.approx((one[0] + two[0], one[1] + two[1]),
                                   rel=1e-12)


def test_bidisperse_crowder_fractions_match_generator_design():
    """
    Generator counts for phi_c 0.25 at pctc 40 give phi_c_one 0.10,
    phi_c_two 0.15, and a type-one share of 40 percent.
    """
    parser = SumRuleCubHeteroRingBiDisperse(
        'al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005bdump2000adump5000',
        'ensemble_long', 'bug')

    fractions = (parser.phi_bulk_c_one, parser.phi_bulk_c_two)
    share = 100 * parser.phi_bulk_c_one / parser.phi_bulk_c

    # phi_i = n_i pi d_i^3 / (6 L^3) with L = 90; the generator rounded each
    # count to an integer, a relative error of at most 0.5/1671 = 3e-4, so
    # at most 4.5e-5 in phi_c_two and 0.012 in the percentage.
    assert fractions == pytest.approx((0.10, 0.15), abs=1e-4)
    assert share == pytest.approx(40.0, abs=0.02)


@pytest.mark.parametrize('artifact, lineage', [
    pytest.param('al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005bdump2000'
                 'adump5000ens1.j20.ring.all', 'segment', id='segment'),
    pytest.param('al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005bdump2000'
                 'adump5000ens1.ring.all', 'whole', id='whole'),
    pytest.param('al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005bdump2000'
                 'adump5000', 'ensemble_long', id='ensemble_long'),
    pytest.param('al4nl5ns400ac1ab5pctc40nc139229nb1671', 'ensemble',
                 id='ensemble'),
    pytest.param('al4nl5ns400ac1ab5pctc40', 'space', id='space'),
])
def test_bidisperse_pct_two_complements_pct_one(artifact, lineage):
    """Every lineage sets pct_two = 100 - pct_one, here 60 for pctc 40."""
    parser = SumRuleCubHeteroRingBiDisperse(artifact, lineage, 'all')

    # 100 - 40.0 is exact in float64.
    assert parser.pct_two == 60.0


@pytest.mark.parametrize('artifact, lineage', [
    pytest.param('al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005bdump2000'
                 'adump5000ens1.j20.ring.all', 'segment', id='segment'),
    pytest.param('al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005bdump2000'
                 'adump5000ens1.ring.all', 'whole', id='whole'),
    pytest.param('al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005bdump2000'
                 'adump5000', 'ensemble_long', id='ensemble_long'),
])
def test_bidisperse_masses_are_size_cubed(artifact, lineage):
    """At a long lineage each computed mass is its size cubed: 4^3, 1, 5^3."""
    parser = SumRuleCubHeteroRingBiDisperse(artifact, lineage, 'all')

    actual = (parser.mmon_large, parser.mcrowd_one, parser.mcrowd_two)

    # Cubes of small whole numbers are exact in float64.
    assert actual == (64.0, 1.0, 125.0)


@pytest.mark.parametrize('artifact, lineage', [
    pytest.param('al4nl5ns400ac1ab5pctc40nc139229nb1671', 'ensemble',
                 id='ensemble'),
    pytest.param('al4nl5ns400ac1ab5pctc40', 'space', id='space'),
])
def test_bidisperse_short_lineages_set_no_masses(artifact, lineage):
    """Ensemble and space set no computed masses, not even a placeholder."""
    parser = SumRuleCubHeteroRingBiDisperse(artifact, lineage, 'all')

    actual = (hasattr(parser, 'mmon_large'), hasattr(parser, 'mcrowd_one'),
              hasattr(parser, 'mcrowd_two'))

    assert actual == (False, False, False)


def test_bidisperse_sets_no_eps_others():
    """The LJ energy unit is documented, not stored as an attribute."""
    parser = SumRuleCubHeteroRingBiDisperse(
        'al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005bdump2000adump5000'
        'ens1.j20.ring.all', 'segment', 'all')

    assert not hasattr(parser, 'eps_others')


@pytest.mark.parametrize(
    'parser_class, group, child_lineage, child_name, parent_lineage',
    [
        pytest.param(TwoMonDepCub, 'bug', 'segment',
                     'am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100'
                     'tdump500ens3.j01.bug', 'whole',
                     id='TwoMonDepCub-segment-whole'),
        pytest.param(TwoMonDepCub, 'bug', 'whole',
                     'am5.0nm2ac1.0nc100hl10.0sd2.0dt0.01bdump50adump100'
                     'tdump500ens3.bug', 'ensemble_long',
                     id='TwoMonDepCub-whole-ensemble_long'),
        pytest.param(SumRuleCyl, 'all', 'segment',
                     'N200epsilon5.0r10.5lz25.0sig2.0nc1000dt0.005bdump2000'
                     'adump5000ens1.j01', 'whole',
                     id='SumRuleCyl-segment-whole'),
        pytest.param(SumRuleCyl, 'all', 'whole',
                     'N200epsilon5.0r10.5lz25.0sig2.0nc1000dt0.005bdump2000'
                     'adump5000ens1', 'ensemble_long',
                     id='SumRuleCyl-whole-ensemble_long'),
        pytest.param(SumRuleCubHeteroRing, 'all', 'whole',
                     'al5nl5ml125ns400ac3nc7735l45dt0.005bdump2000adump5000'
                     'ens4.ring.all', 'ensemble_long',
                     id='SumRuleCubHeteroRing-whole-ensemble_long'),
        pytest.param(SumRuleCubHeteroLinear, 'all', 'whole',
                     'al6nl5ml216ns400ac3nc7735l45dt0.005bdump2000adump5000'
                     'ens1.linear.all', 'ensemble_long',
                     id='SumRuleCubHeteroLinear-whole-ensemble_long'),
        pytest.param(TransFociCyl, 'bug', 'segment',
                     'epss5epsl5r10.5al3nl5ml27ns400ac1nc8340lz69.5dt0.005'
                     'bdump2000adump5000ens8.j05.ring.bug', 'whole',
                     id='TransFociCyl-segment-whole'),
        pytest.param(TransFociCub, 'bug', 'whole',
                     'al5.0nl5ml1.0ns400ac1.0nc100l15.0dt0.005bdump2000'
                     'adump5000ens1', 'ensemble_long',
                     id='TransFociCub-whole-ensemble_long'),
        pytest.param(HnsCub, 'nucleoid', 'segment',
                     'N200kbmm2nh12ac2l25epshc1nc11937ens2.j5.ring.nucleoid',
                     'whole', id='HnsCub-segment-whole'),
        pytest.param(HnsCyl, 'all', 'segment',
                     'N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152ens1.j01.ring.all',
                     'whole', id='HnsCyl-segment-whole'),
        pytest.param(HnsCyl, 'all', 'whole',
                     'N200kbmm2r4.5nh4ac1lz75epshc1.0nc1152ens1.ring.all',
                     'ensemble_long', id='HnsCyl-whole-ensemble_long'),
        pytest.param(SumRuleCubHeteroRingBiDisperse, 'all', 'segment',
                     'al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005'
                     'bdump2000adump5000ens1.j20.ring.all', 'whole',
                     id='SumRuleCubHeteroRingBiDisperse-segment-whole'),
        pytest.param(SumRuleCubHeteroRingBiDisperse, 'bug', 'whole',
                     'al4nl5ns400ac1ab5pctc40nc139229nb1671l45dt0.005'
                     'bdump2000adump5000ens1.ring.bug', 'ensemble_long',
                     id='SumRuleCubHeteroRingBiDisperse-whole-ensemble_long'),
    ],
)
def test_parser_rebuilt_parent_name_keeps_attributes(
        parser_class, group, child_lineage, child_name, parent_lineage):
    """A rebuilt parent name parses to its child's attribute values."""
    child = parser_class(child_name, child_lineage, group)
    parent = parser_class(getattr(child, parent_lineage), parent_lineage,
                          group)
    values = attrgetter(*parser_class._genealogy_attributes[parent_lineage])

    # Both sides parse the same digits, so they agree to float64 round-off.
    assert values(parent) == pytest.approx(values(child), rel=1e-12,
                                           abs=1e-12)


@pytest.mark.parametrize(
    'parser_class, group, lineage, artifact, parent_lineage, expected',
    [
        pytest.param(SumRuleCyl, 'all', 'segment',
                     'N200epsilon5.0r10.5lz25.0sig2.0nc1000dt0.005bdump2000'
                     'adump5000ens1.j01', 'segment',
                     'N200epsilon5.0r10.5lz25.0sig2.0nc1000dt0.005bdump2000'
                     'adump5000ens1.j01', id='SumRuleCyl-segment'),
        pytest.param(SumRuleCyl, 'all', 'segment',
                     'N200epsilon5.0r10.5lz25.0sig2.0nc1000dt0.005bdump2000'
                     'adump5000ens1.j01', 'whole',
                     'N200epsilon5.0r10.5lz25.0sig2.0nc1000dt0.005bdump2000'
                     'adump5000ens1', id='SumRuleCyl-whole'),
        pytest.param(SumRuleCyl, 'all', 'segment',
                     'N200epsilon5.0r10.5lz25.0sig2.0nc1000dt0.005bdump2000'
                     'adump5000ens1.j01', 'ensemble_long',
                     'N200epsilon5.0r10.5lz25.0sig2.0nc1000dt0.005bdump2000'
                     'adump5000', id='SumRuleCyl-ensemble_long'),
        pytest.param(TransFociCyl, 'bug', 'segment',
                     'epss5epsl5r10.5al3nl5ml27ns400ac1nc8340lz69.5dt0.005'
                     'bdump2000adump5000ens8.j05.ring.bug', 'segment',
                     'epss5epsl5r10.5al3nl5ml27ns400ac1nc8340lz69.5dt0.005'
                     'bdump2000adump5000ens8.j05', id='TransFociCyl-segment'),
        pytest.param(TransFociCyl, 'bug', 'whole',
                     'epss5epsl5r10.5al3nl5ml27ns400ac1nc8340lz69.5dt0.005'
                     'bdump2000adump5000ens8.ring.bug', 'whole',
                     'epss5epsl5r10.5al3nl5ml27ns400ac1nc8340lz69.5dt0.005'
                     'bdump2000adump5000ens8', id='TransFociCyl-whole'),
        pytest.param(HnsCub, 'nucleoid', 'segment',
                     'N200kbmm2nh12ac2l25epshc1nc11937ens2.j5.ring.nucleoid',
                     'ensemble_long', 'N200kbmm2nh12ac2l25epshc1nc11937',
                     id='HnsCub-ensemble_long'),
    ],
)
def test_parser_rebuilt_parent_name_is_name_prefix(
        parser_class, group, lineage, artifact, parent_lineage, expected):
    """A long-form parent name is the artifact name cut after its last key."""
    parser = parser_class(artifact, lineage, group)

    assert getattr(parser, parent_lineage) == expected


# Declarations


@pytest.mark.parametrize('parser_class, lineage, expected', [
    pytest.param(SumRuleCyl, 'whole',
                 OrderedDict([('nmon', 'N'), ('epsilon', 'epsilon'),
                              ('dcyl', 'r'), ('lcyl', 'lz'),
                              ('dcrowd', 'sig'), ('ncrowd', 'nc'),
                              ('dt', 'dt'), ('bdump', 'bdump'),
                              ('adump', 'adump'), ('ensemble_id', 'ens')]),
                 id='SumRuleCyl-whole'),
    pytest.param(SumRuleCyl, 'ensemble',
                 OrderedDict([('nmon', 'N'), ('dcyl', 'D'),
                              ('dcrowd', 'ac'), ('ncrowd', 'nc')]),
                 id='SumRuleCyl-ensemble'),
    pytest.param(HnsCub, 'space',
                 OrderedDict([('nmon', 'N'), ('bend_mm', 'kbmm'),
                              ('nhns', 'nh'), ('dcrowd', 'ac'),
                              ('eps_hc', 'epshc')]),
                 id='HnsCub-space'),
])
def test_genealogy_attributes_follow_pattern_order(parser_class, lineage,
                                                   expected):
    """Each lineage maps attributes to keywords in the pattern's order."""
    table = parser_class._genealogy_attributes[lineage]

    assert (type(table), list(table.items())) == (OrderedDict,
                                                  list(expected.items()))


@pytest.mark.parametrize('parser_class', [
    pytest.param(ParserBase, id='ParserBase'),
    pytest.param(type('PatternlessParser', (ParserBase,), {}),
                 id='subclass-without-patterns'),
])
def test_parser_without_patterns_cannot_be_instantiated(parser_class):
    """A parser class that declares no patterns raises TypeError."""
    with pytest.raises(TypeError, match='declares no lineage patterns'):
        parser_class('N200', 'space', 'bug')


def test_declared_parser_parses_its_names(chain_declaration):
    """A class built from a valid declaration parses names of its grammar."""
    chain_cub = type('ChainCub', (ParserBase,), chain_declaration)

    parser = chain_cub('N100ac2.0nc50l10.0ens1.bug', 'whole', 'bug')

    assert (parser.nmon, parser.dcrowd, parser.ncrowd, parser.lcube,
            parser.ensemble_long, parser.ensemble, parser.space) == (
        100, 2.0, 50, 20.0, 'N100ac2.0nc50l10.0', 'N100ac2nc50',
        'N100ac2')


@pytest.mark.parametrize('override, match', [
    pytest.param({'_geometry': None}, 'must declare', id='no-geometry'),
    pytest.param({'_dump_keywords': None}, "'_dump_keywords'",
                 id='no-dump-keywords'),
    pytest.param({'_dump_keywords': {'bug': 'bdump'}},
                 "_dump_keywords must have one entry per group",
                 id='group-no-dump-keyword'),
    pytest.param(
        {'_dump_keywords': {'bug': 'bdump', 'all': 'adump',
                            'nucleoid': 'ndump'}},
        "_dump_keywords must have one entry per group",
        id='dump-keyword-of-undeclared-group'),
    pytest.param(
        {'_dump_keywords': {'bug': 'tdump', 'all': 'adump'}},
        "the dump keyword 'tdump' of group 'bug' is not set by the "
        "'segment' pattern or '_defaults'",
        id='dump-keyword-not-set'),
    pytest.param(
        {'_patterns': {'segment': ('N', 'ac', 'nc', 'l', 'ens', 'j'),
                       'whole': ('N', 'ac', 'nc', 'l', 'ens'),
                       'ensemble_long': ('N', 'ac', 'nc', 'l'),
                       'ensemble': ('N', 'ac', 'nc')}},
        "_patterns must have one entry per lineage", id='lineage-no-pattern'),
    pytest.param(
        {'_project_attributes': {'segment': [], 'whole': []}},
        "_project_attributes must have one entry per lineage",
        id='lineage-no-project-attributes'),
    pytest.param(
        {'_patterns': {'segment': ('N', 'ac', 'nc', 'l', 'ens', 'j'),
                       'whole': ('N', 'ac', 'nc', 'l', 'ens'),
                       'ensemble_long': ('N', 'ac', 'nc', 'l'),
                       'ensemble': ('N', 'ac', 'nc'),
                       'space': ('N', 'xx')}},
        "keyword 'xx' in the 'space' pattern has no KeywordSpec",
        id='undeclared-keyword'),
    pytest.param(
        {'_patterns': {'segment': ('N', 'ac', 'nc', 'l', 'ens', 'j'),
                       'whole': ('N', 'ac', 'nc', 'l', 'ens'),
                       'ensemble_long': ('N', 'ac', 'nc', 'l'),
                       'ensemble': ('N', 'ac', 'nc'),
                       'space': ('N', 'ac', 'N')}},
        "attribute 'nmon' appears more than once in the 'space' pattern",
        id='duplicate-attribute'),
    pytest.param(
        {'_totals': {'all': ('m', 'x')}},
        "total 'all' sums 'x', which is not a declared species",
        id='total-of-undeclared-species'),
    pytest.param(
        {'_patterns': {'segment': ('N', 'ac', 'nc', 'l', 'ens', 'j'),
                       'whole': ('N', 'ac', 'nc', 'l', 'ens'),
                       'ensemble_long': ('N', 'ac', 'nc'),
                       'ensemble': ('N', 'ac', 'nc'),
                       'space': ('N', 'ac')}},
        "the bulk values need 'lcube', which the 'ensemble_long' pattern",
        id='bulk-input-not-parsed'),
])
def test_invalid_declaration_raises_type_error_at_class_creation(
        chain_declaration, override, match):
    """An incomplete or inconsistent declaration fails when it is defined."""
    with pytest.raises(TypeError, match=re.escape(match)):
        type('BadParser', (ParserBase,), {**chain_declaration, **override})


# KeywordSpec


@pytest.mark.parametrize('spec, token, expected', [
    pytest.param(KeywordSpec('nmon', int), '200', 200, id='int'),
    pytest.param(KeywordSpec('ncrowd', int), '1000.', 1000,
                 id='int-trailing-dot'),
    pytest.param(KeywordSpec('dcrowd', float), '2.5', 2.5, id='float'),
    pytest.param(KeywordSpec('lcube', float, scale=2.0), '25.0', 50.0,
                 id='half-side'),
    pytest.param(KeywordSpec('dcyl', float, scale=2.0, offset=-1.0), '10.5',
                 20.0, id='wall-radius'),
])
def test_keyword_spec_converts_token(spec, token, expected):
    """A token converts to `scale * value + offset`, typed as declared."""
    actual = spec.parse(token)

    # The conversions are exact in float64 for these tokens.
    assert (actual, type(actual)) == (pytest.approx(expected, rel=1e-12),
                                      type(expected))


@pytest.mark.parametrize('spec, token, match', [
    pytest.param(KeywordSpec('dcrowd', float), '.', "'.' is not a number",
                 id='dot'),
    pytest.param(KeywordSpec('dcrowd', float), '', "'' is not a number",
                 id='empty'),
    pytest.param(KeywordSpec('nmon', int), '1.5', "'1.5' is not an integer",
                 id='non-integral'),
])
def test_keyword_spec_rejects_bad_token(spec, token, match):
    """A token that is not a number, or not whole for int, is rejected."""
    with pytest.raises(ValueError, match=re.escape(match)):
        spec.parse(token)


@pytest.mark.parametrize('kwargs, match', [
    pytest.param({'attr': 'name', 'kind': str}, "must be int or float",
                 id='str-kind'),
    pytest.param({'attr': 'nmon', 'kind': int, 'offset': 1.0},
                 "only a float keyword", id='int-with-offset'),
])
def test_keyword_spec_rejects_bad_declaration(kwargs, match):
    """A KeywordSpec of another type, or an int with a conversion, fails."""
    with pytest.raises(TypeError, match=match):
        KeywordSpec(**kwargs)
