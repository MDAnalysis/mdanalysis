import pooch

from MDAnalysis.core.universe import Universe


def fetch_pdb(PDB_IDS=None,
            cache_path=None,
            progress_bar=False,
            file_format="pdb.gz",
            ):
    
    ## Have to do this approach instead of Pooch.retrieve in order to prevent known hash warning from showing up
    
    # Handles the case if input is a string
    if isinstance(PDB_IDS, str):
        PDB_IDS = (PDB_IDS,)

    # Handles multiple tuples and lists
    registry_dictionary = {f'{PDB_ID}.{file_format}': None for PDB_ID in PDB_IDS}
    filepath_list = []
    universe_list = []
    
    downloader = pooch.create(
        path=cache_path,
        base_url="https://files.rcsb.org/download/",
        registry=registry_dictionary
    )

    # Make so that in the case of passing one string- the output is just one universe
    for file_name in registry_dictionary.keys():
        filepath_list.append(downloader.fetch(fname=file_name, progressbar=progress_bar))

    for path in filepath_list:
        universe_list.append(Universe(path))
        
    return universe_list





u = fetch_pdb(PDB_IDS='1AKE', cache_path='~/tmp', progress_bar=True)

print(u)


print(fetch_pdb(PDB_IDS=['1AKE', '4BWZ'] , cache_path='~/tmp', progress_bar=True))
