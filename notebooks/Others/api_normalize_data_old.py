

    
    
def rescale_by_phi_c_0(
    group: pd.DataFrame,
    prop: str,
    space: str
) -> pd.Series:
    """
    Rescale a given property within a group based on its value where
    'phi_c_bulk_round' is zero.

    Parameters
    ----------
    group : pd.DataFrame
        A subset of the dataframe grouped by a particular column
        (e.g., 'space').
    prop : str
        The name of the property column which needs to be rescaled.
    space : str
        The space type, used for warnings.

    Returns
    -------
    pd.Series
        A series representing the rescaled property values for the group.

    Notes
    -----
    If the value of the property at 'phi_c_bulk_round' equals 0 or NaN, a
    warning is issued and NaN values are returned for the entire group
    """
    # Find the value of the property where phi_c_bulk_round is 0
    print(group['space'].iloc[0])
    value_at_conditions = \
        group[group['phi_c_bulk_round'] == 0][f"{prop}-mean"].iloc[0]
    print(value_at_conditions)
    # If the value is 0 or NaN, avoid division by zero
    if value_at_conditions == 0 or pd.isna(value_at_conditions):
        warnings.warn(
            f"The '{prop}' value in the absence of crowders "
            f"(phi_c=0) for space '{space}', the values of '{prop}-norm' "
            "at all the values of phi_c are set to 'np.nan'."
        )
        return pd.Series([np.nan] * len(group))
    return group[f"{prop}-mean"] / value_at_conditions


def normalize_data(
    df: pd.DataFrame,
    project: str,
    properties: List[str]
) -> pd.DataFrame:
    """
    Normalize the data in the dataframe by specific conditions.

    For the given properties in the dataframe, this function calculates their
    normalized values based on their values where 'phi_c_bulk_round' is zero.

    Parameters
    ----------
    df : pd.DataFrame
        The input dataframe containing the data.
    project : str
        The name of the project, used for additional normalizations.
    properties : List[str]
        A list of property columns which should be normalized.

    Returns
    -------
    pd.DataFrame
        A dataframe with additional columns representing normalized values.

    Notes
    -----
    Additional normalization is performed for specific projects:
        ['HnsCyl', 'HnsCub'].
    """
    df_copy = df.copy()
    norm_props = [
        prop.split('-')[0] for prop in properties if prop.endswith('mean')]
    for prop in norm_props:
        series = df_copy.groupby('space').apply(
            lambda group: rescale_by_phi_c_0(
                group, prop, group['space'].iloc[0])).reset_index(
                    level=0, drop=True)
        print('series:', series)
        df_copy[f"{prop}-norm"] = series

    if project in ['HnsCyl', 'HnsCub']:
        n_patch_per_cor = 2
        hpatch_cols = ['nBoundHnsPatch', 'nFreeHnsPatch', 'nEngagedHnsPatch']
        hcore_cols = ['nFreeHnsCore', 'nBridgeHnsCore', 'nDangleHnsCore',
                      'nCisHnsCore', 'nTransHnsCore']
        for col in hpatch_cols:
            df_copy[col + '-norm'] = \
                df_copy[col + '-mean'] / (df_copy['nhns'] * n_patch_per_cor)
        for col in hcore_cols:
            df_copy[col + '-norm'] = df_copy[col + '-mean'] / df_copy['nhns']

    return df_copy