import pandas as pd
import click
import numpy as np

# @click.command()
# @click.argument('df')
# @click.argument('taxa_column')
# @click.argument('delimeter')
# @click.argument('taxa_level')
# @click.option('--output_fp', '-o', default = 'genome.bed')

def summarize_taxa(df, count_header, taxa_header, taxa_column, delimeter, taxa_level, output_fp=None):
    
    """Summarize a long-form count table at a given taxanomic level. This code needs some more work to make it 
       more useful and generic (not dependent on column headers)
    
    Args:
        df (pandas dataframe)         : input data frame with columns indicating samples id (sid)
        count_header (string)         : name of column with counts to sum
        taxa_header (string)          : name of column with the taxa to parse
        taxa_column (list of strings) : list of the taxa levels
        delimeter (character)         : character that deliminates taxa levels   
        taxa_level (int)              : level to summarize
        output_fp (str)               : output file path (optional)
        
    
    Yields:
        sum_taxa: data frame where counts from the same level and sample are summed 
    
    """
    
    new_level = taxa_column[0:taxa_level]

    # Not reproducibleundil the column names are 
    df[taxa_column] = df.OTU.str.split(delimeter, expand=True)

    # Add a true/false arguement here
    df = df[df['domain']!="Eukaryota"]

    new_level.insert(0,'sid')

    taxa_sum = df.groupby(new_level)[count_header].sum().reset_index()
    taxa_sum[taxa_header] = taxa_sum[taxa_column[0:taxa_level]].apply(lambda row: '|'.join(row.values.astype(str)), axis=1)
    
    return taxa_sum
    
    if __name__ == '__main__':
        summarize_taxa()