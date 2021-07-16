#!/usr/bin/python3
import re



with open('bp-sh(1).tags', 'r') as required_rsid:
    list_rsid = [line.strip() for line in required_rsid]
  
    
with open('1kg_503_eur_82M_ref.bim', mode='r') as in_file:
    for line in in_file:
        list_line = re.findall(r'\w+', line)
        if list_line[1] in list_rsid:
            with open('required_rsid.bim', mode='a') as out_file:
                new_line = list_line[1] + ',' + list_line[0] + ',' + list_line[3] + ',' + list_line[4] + ',' + list_line[5]
                out_file.write(f'{new_line}\n')
