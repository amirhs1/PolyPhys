import datetime

from ..polyphys import CellAttributes

def data_file_generator(template, fnames):
    
    if not template[0].endswith('data'):
        raise Exception("The template file is not Lammps data file.")
    
    for fname in fnames:
        if not fname.endswith('lammpstrj'):
            raise Exception("The input files are not Lammps trajectroy files.")
    
    today = datetime.date.today().strftime('%Y%m%d')
    firstline = f"LAMMPS data file via PipeLine.data_file_generator, based on LAMMPS (template) data file via write_data, version {today}, timestep = 71000000"
    
    for fname in fnames:
        output = fname.split('/')[-1].split('.lammpstrj')[0]+'.data'
        cell_attrs = CellAttributes(fname,'cylindrical')
        xmax = (cell_attrs.dcyl+1.0) / 2.0
        ymax = (cell_attrs.dcyl+1.0) / 2.0
        zmax = cell_attrs.zlen / 2.0
        with open(template[0],'r') as data_in,\
        open(output,'w') as data_out:
            _ = data_in.readline()
            _ = data_in.readline()
            data_out.write(firstline)
            data_out.write('\n')
            data_out.write('\n')
            line = data_in.readline()
            data_out.write(line)
            line = data_in.readline()
            data_out.write(line)
            line = data_in.readline()
            data_out.write(line)
            line = data_in.readline()
            data_out.write(line)
            _ = data_in.readline()
            _ = data_in.readline()
            _ = data_in.readline()
            _ = data_in.readline()
            data_out.write('\n')
            data_out.write(str(-1.0*xmax)+' '+str(xmax)+' xlo xhi \n')
            data_out.write(str(-1.0*ymax)+' '+str(ymax)+' ylo yhi \n')
            data_out.write(str(-1.0*zmax)+' '+str(zmax)+' zlo zhi \n')
            line = data_in.readline()
            while(line):
                data_out.write(line)
                line = data_in.readline()