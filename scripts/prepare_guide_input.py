import os, sys
from pathlib import Path
import pandas as pd

def prepare_guide_input(guide_file):
    dir_path = os.path.dirname(guide_file)
    guide_file_name = os.path.basename(guide_file).removesuffix('.xlsx').removesuffix('.xls')
    output_file = os.path.join(dir_path, f"{guide_file_name}.csv")
    df = pd.read_excel(guide_file, sheet_name=0, header=0)
    df_sub = df[['ID', 'Spacer']]
    df_sub.columns = ['guide_id', 'guide sequence']
    df_sub.to_csv(output_file, sep="\t", index=False)

    


def main():
    guide_file = sys.argv[1]
    prepare_guide_input(guide_file)


if __name__ == "__main__":
    main()