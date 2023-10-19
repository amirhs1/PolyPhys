from glob import glob
from polyphys.manage.parser import HnsCub
from polyphys.probe.prober import stamps_report

stamps = glob("./N*/N*/N*stamps.csv")
for stp in stamps:
    sim_info = HnsCub(
        stp,
        'whole',
        'cubic',
        'nucleoid',
        'ring'
    )
    sim_name = sim_info.lineage_name + "-" + sim_info.group
    n_frames = 20001

    # Simulation stamps:
    outfile = stp.split("-stamps.csv")[0]
    outfile = outfile + "-stamps-new.csv"
    stamps_report(outfile, sim_info, n_frames)
