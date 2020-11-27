import os
import sys
import re
sys.path.append('../..')

re_pattern1 = r'# Dictionary of parameter names'
rg_pattern1 = re.compile(re_pattern1,re.ASCII)

str1 = "# List of parameter dictionaries with names, initial values,"
str2 = "# and min/max bounds. Set 'vary': False to hold a param constant."
str3 = "    :param kws: dict of optional args (eg load, fsf, zsf)"
str4 = "def model(w, params, **kws):"

file_list = os.listdir()
models = []
for file in file_list:
    if file.endswith(".py") and not ("ModelConvert" in file):
        models.append(file)

new_lines = []
for file in models:
    print(file)
    revised = []
    with open(file, mode='r', encoding='utf-8', newline=None) as f:
        for line in f:
            g = rg_pattern1.search(line)
            if g:
                revised.append(str1 + '\n')
                revised.append(str2 + '\n')
                next(f)
                next(f)
                next(f)
                continue
            if 'model(w, load, params)' in line:
                line = str4 + '\n'
            if ':param load:' in line:
                line = next(f)
            if ':param params:' in line:
                add_line = str3 + '\n'
            else:
                add_line = None
            revised.append(line)
            if add_line:
                revised.append(add_line)

    with open(file, mode='w', encoding='utf-8', newline=None) as f:
        for line in revised:
            f.write(line)
