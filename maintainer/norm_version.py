from packaging.version import Version


def norm_version(version_str: str):
    """
    Normalize an input version string in the same way that setuptools' dist
    does.

    Parameters
    ----------
    version_str : str
      A version string to normalize.

    Returns
    -------
    str
      Normalised version string
    """
    return str(Version(version_str))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--file', type=str, help="file with version to parse")
    args = parser.parse_args()

    with open(args.file) as filed:
        ver = filed.readlines()[0].strip("\n")

    print(norm_version(version_str=ver))
