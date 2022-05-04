import os, glob


def read_files(directory: str) -> dict[str, list[str]]:
    '''
    Read the files in the given directory and return how many and which
    platforms

    Parameters
    ----------
    directory: str
        Directory where the files are located

    Returns
    ----------
    platform_files: dict[str, list[str]]
        Types of platforms found and their files associated
    '''
    os.chdir(directory)
    files = glob.glob('*')
    p180_files = []
    nmr_files = []
    for f in files:
        if 'ADMCDUKEP180' in f:
            p180_files.append(f)
        elif 'ADNINIGHTINGALE2' in f:
            nmr_files.append(f)
    platform_files = {'p180': p180_files,
                      'nmr': nmr_files}
    return(platform_files)
