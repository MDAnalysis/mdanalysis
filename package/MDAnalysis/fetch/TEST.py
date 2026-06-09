#!/usr/bin/env python3
## File for testing changes 
## Will get deleted when PR is DONE


# import MDAnalysis as mda

# from MDAnalysis.fetch.fetchers import StaticFetcher
# from MDAnalysis.fetch.pdb import from_PDB

# import hashlib


# print(hashlib.algorithms_available)
# print(" ")
# print(hashlib.algorithms_guaranteed)


# s = StaticFetcher(reuse_connection=False, hash='sha256')

# print(s.cache_path)


from MDAnalysis.fetch.fetchers import StaticFetcher

s = StaticFetcher()

s.fetch()