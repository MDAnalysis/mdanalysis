import pooch

from MDAnalysis.core.universe import Universe


def fetch_pdb(PDB_IDS=None,
            cache_path=None,
            progressbar=False,
            file_format="pdb.gz",
            ):
    # note, progress_bar required tdqm. make an optional depedency if not already in
    
    ## Have to do this approach instead of Pooch.retrieve in order to prevent known hash warning from showing up

    # Handles the case if input arguement is a string (i.e. only one pdb)
    if isinstance(PDB_IDS, str):
        PDB_IDS = (PDB_IDS,)

    # Handles multiple tuples and lists
    registry_dictionary = {f'{PDB_ID}.{file_format}': None for PDB_ID in PDB_IDS}
    
    
    downloader = pooch.create(
        path=cache_path,
        base_url="https://files.rcsb.org/download/",
        registry=registry_dictionary
    )

<<<<<<< HEAD
    if len(PDB_IDS) == 1:
        return downloader.fetch(fname=tuple(registry_dictionary.keys())[0], progressbar=progressbar)
=======
    print(registry_dictionary.items())


    if len(PDB_IDS) == 1: # Ensures only one Path is returned in the case of one PDB_ID input argument
        for file_name in registry_dictionary.keys():
            return downloader.fetch(fname=file_name, progressbar=progressbar)
>>>>>>> parent of 03638c8ef (Cleaned up return logic with syntactic sugar)
    else:
        filepath_list = []
        for file_name in registry_dictionary.keys():
            filepath_list.append(downloader.fetch(fname=file_name, progressbar=progressbar))
        return filepath_list



print(fetch_pdb(PDB_IDS='1AKE', cache_path='~/tmp', progressbar=True))
print(fetch_pdb(PDB_IDS=['1AKE', '4BWZ'] , cache_path='~/tmp', progressbar=True))
