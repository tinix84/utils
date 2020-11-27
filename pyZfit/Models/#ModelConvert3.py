import os
import sys
import re
sys.path.append('../..')

# General re for time only
re_pattern = r'(.*params\[.*\])\.value'         # 'params' line without '.value' in group 1
rg_pattern = re.compile(re_pattern,re.ASCII)

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
            g = rg_pattern.search(line)
            if g:
                new_lines.append(g.group(1) + '\n')
                revised.append(g.group(1) + '\n')
            else:
                revised.append(line)

    with open(file, mode='w', encoding='utf-8', newline=None) as f:
        for line in revised:
            f.write(line)

with open("new_lines.txt", mode="w") as f:
    for line in new_lines:
        f.write(line)