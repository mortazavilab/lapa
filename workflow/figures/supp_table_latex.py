import pandas as pd


for i, table in enumerate(snakemake.input['tables']):
    df = pd.read_csv(table)

    with open(f'reports/tables/supp_table{i + 1}.txt', 'w') as f:
        f.write(df.to_latex(index=False))
