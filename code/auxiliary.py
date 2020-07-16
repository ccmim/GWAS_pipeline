
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
