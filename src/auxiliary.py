import os
import yaml

def is_yaml_file(x):
    if isinstance(x, str):
        return x.endswith("yaml") or x.endswith("yml")
    return False

# copied from http://thoughtsbyclayg.blogspot.com/2008/10/parsing-list-of-numbers-in-python.html
# with a few modifications
def parseIntSet(nputstr=""):

    '''
    :param nputstr: set of comma-separated number, number ranges (such as 1-22) or X,Y
    :return:
    '''

    selection = set()
    invalid = set()
    # tokens are comma separated values

    tokens = [x.strip() for x in nputstr.split(',')]
    for i in tokens:
        try:
            # autosomal chromosomes are integers between 1 and 22
            if int(i) <= 22 and int(i) >= 1:
                selection.add(i)
            else:
                invalid.add(i)
        except:
            # if the token is not a number, then it might be a range, like 5-9
            try:
                token = [int(k.strip()) for k in i.split('-')]
                if len(token) > 1:
                    token.sort()
                    # we have items separated by a dash
                    # try to build a valid range
                    first = token[0]
                    last = token[len(token) - 1]
                    for x in range(first, last + 1):
                        selection.add(str(x))
            except:
                if i == "X" or i == "Y":
                    selection.add(i)
                else:
                    # if not an int nor a range nor X/Y...
                    invalid.add(i)

    # Report invalid tokens before returning valid selection
    # print("Invalid set: " + str(invalid))

    return selection

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def get_repo_rootdir():
    import shlex
    from subprocess import check_output
    repo_rootdir = check_output(shlex.split("git rev-parse --show-toplevel")).strip().decode('ascii')
    return repo_rootdir


def unfold_config(token, no_unfolding_for=[]):
    '''
    Parameters: 
      token: a recursive structure composed of a path to a yaml file or a dictionary composed of such structures.
      no_unfolding_for: a list of dict keys for which the yaml shouldn't be unfolded, and instead kept as a path
    Returns: A dictionary with all the yaml files replaces by their content.
    '''
    repo_rootdir = get_repo_rootdir()
    yaml_dir = os.path.join(repo_rootdir, "config_files")
    if is_yaml_file(token):
        #TODO: COMMENT AND DOCUMENT THIS!!!
        try:
            token = yaml.safe_load(open(token))
        except FileNotFoundError:
            kk = open(os.path.join(yaml_dir, token))
            token = yaml.safe_load(kk)
    if isinstance(token, dict):
        for k, v in token.items():
            if k not in no_unfolding_for:
                token[k] = unfold_config(v, no_unfolding_for)
    return token