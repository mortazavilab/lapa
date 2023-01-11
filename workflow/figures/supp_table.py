import pandas as pd


writer = pd.ExcelWriter(snakemake.output['excel'])

for i, table in enumerate(snakemake.input['tables']):
    df = pd.read_csv(table)
    df.to_excel(writer, index=False,
                sheet_name=f'supplementary table {i + 1}')

writer.save()
