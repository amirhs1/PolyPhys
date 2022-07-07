from glob import glob
from polyphys.construct.builder import bug_data_file_generator


def main() -> None:
    path = "./*.bug.lammpstrj"
    bug_trjs = glob(path)
    data_template = glob("./data_template-biaxial-bug-N*.data")
    data_template = data_template[0]
    bug_data_file_generator(
        data_template,
        bug_trjs,
        geometry='biaxial',
        lineage='whole'
    )


main()
