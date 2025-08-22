import pooch

from MDAnalysis.core.universe import Universe


def fetch_pdb(PDB_IDS=None,
            cache_path=None,
            progressbar=False,
            file_format="pdb.gz",
            ):
    
    # Have to do this dictionary approach instead of using Pooch.retrieve in order to prevent the hardcoded known hash warning from showing up
    
    if isinstance(PDB_IDS, str):
        PDB_IDS = (PDB_IDS,)

    registry_dictionary = {f'{PDB_ID}.{file_format}': None for PDB_ID in PDB_IDS}
  
    downloader = pooch.create(
        path=cache_path,
        base_url="https://files.rcsb.org/download/",
        registry=registry_dictionary
    )

    if len(PDB_IDS) == 1:
        return str(downloader.fetch(fname=tuple(registry_dictionary.keys())[0], progressbar=progressbar))
    else:
        return [downloader.fetch(fname=file_name, progressbar=progressbar) for file_name in registry_dictionary.keys()]


print(fetch_pdb(PDB_IDS='1AKE', cache_path='~/tmp', progressbar=True))
print(fetch_pdb(PDB_IDS=['1AKE', '4BWZ'] , cache_path='~/tmp', progressbar=True))
