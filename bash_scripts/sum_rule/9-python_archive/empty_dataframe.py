def empty_dataframe(df_csv):
    """
    empty_dataframe create an empty distributionsin which all the frequencies are zero based on the csv file of an input distribution.
    
    Parameters:
    df_csv (str): address of the dataframe
    
    Return:
    an empty distribtuions (pandas dataframe) in which all the frequencies are zero.
    
    Requirements:
    Pandas
    """
    df = pd.read_csv(df_csv,header=0,index_col=0)
    for col in df.columns:
        df[col].values[:] = 0.0
    return df