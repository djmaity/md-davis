import os

def dir_path(path_string):
    if os.path.isdir(path_string):
        return path_string
    else:
        raise NotADirectoryError(path_string)


def file_path(path_string):
    if os.path.isfile(path_string):
        return path_string
    else:
        raise FileNotFoundError(path_string)
