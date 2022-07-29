# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.7
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
from bs4 import BeautifulSoup
import numpy as np
import os
import logging
import json
import argparse


# %%
def format_key_str(s):
    out = s
    l = ['#','\'',' ','(',')','>=','%']
    for i in l:
        out = out.replace(i,'')
    return out


def get_total_report(html_file):
    outdict = {}
    with open(html_file,'rt') as h:
        soup = BeautifulSoup(h, 'html.parser')
    tmpout = soup.find_all(id='total-report-json')[0].string
    tmpout = json.loads(tmpout.strip())
#     print(tmpout.keys())
    for info in tmpout['report']:
        if info[0] == 'Mismatches':
            for mis in info[1]:
                tmp_key = format_key_str(mis['metricName'])
                outdict[tmp_key] = mis['values'][0]
 
    for info in tmpout['report']:
        if info[0] in ['Mismatches', 'Statistics without reference']:
            for mis in info[1]:
                tmp_key = format_key_str(mis['metricName'])
                outdict[tmp_key] = mis['values'][0] 
    return outdict



# %%
def get_contigs_lengths(html_file):
    outdict = {}
    with open(html_file,'rt') as h:
        soup = BeautifulSoup(h, 'html.parser')
    tmpout = soup.find_all(id='contigs-lengths-json')[0].string
    tmpout = json.loads(tmpout.strip())
#     print(tmpout.keys())
#     outdict['lists_of_lengths'] = tmpout['lists_of_lengths'][0]
    outdict = list(np.cumsum(tmpout['lists_of_lengths'][0]))
    return outdict



# %%

# %%
def get_coord_nx(html_file):
    outdict = {}
    with open(html_file,'rt') as h:
        soup = BeautifulSoup(h, 'html.parser')
    tmpout = soup.find_all(id='coord-nx-json')[0].string
    tmpout = json.loads(tmpout.strip())
#     print(tmpout.keys())
    outdict['coord_y'] = tmpout['coord_y'][0]
    outdict['coord_x'] = tmpout['coord_x'][0]
    return outdict



# %%

# %%
def get_gc(html_file):
    outdict = {}
    with open(html_file,'rt') as h:
        soup = BeautifulSoup(h, 'html.parser')
    tmpout = soup.find_all(id='gc-json')[0].string
    tmpout = json.loads(tmpout.strip())
#     print(tmpout.keys())
    outdict['list_of_GC_distributions'] = {}
    outdict['list_of_GC_distributions']['coord_x'] = tmpout['list_of_GC_distributions'][0][0]
    outdict['list_of_GC_distributions']['coord_y'] = tmpout['list_of_GC_distributions'][0][1]
    
    outdict['list_of_GC_contigs_distributions'] = {}
    outdict['list_of_GC_contigs_distributions']['coord_x'] = tmpout['list_of_GC_contigs_distributions'][0][0]
    outdict['list_of_GC_contigs_distributions']['coord_y'] = tmpout['list_of_GC_contigs_distributions'][0][1]
    
    return outdict



# %%

# %%

# %%
def parse_quast_html(html_file):
    outdict = {}
    outdict['sta'] = get_total_report(html_file)
    outdict['contig_length'] = get_contigs_lengths(html_file)
    outdict['Nx'] = get_coord_nx(html_file)
    outdict['gc_content'] = get_gc(html_file)
    return outdict


# %%

# %%
if __name__ == '__main__':
    parse = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parse.add_argument('-i', '--input', required=True, help='quast report output html file for format')
    parse.add_argument('-o', '--output', default=None, help='output json, default stdout')
    args = parse.parse_args()
    
    logging.basicConfig(level='INFO')
    html_file = args.input
    outdict = parse_quast_html(html_file)
    if args.output:
        try:
            with open(args.output,'w') as h:
                json.dump(outdict, h, indent=2)
        except Exception as e:
            logging.error(e)
    else:
        print(outdict)
