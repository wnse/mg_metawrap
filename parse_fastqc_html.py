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
import base64
import pandas as pd
import re
import base64
import os
import logging
import json
import argparse


# %%

# %%
def get_fqc(html_file):
    outdict = {}
    outdict['html'] = os.path.realpath(html_file)
    outdict['sta'] = {}
    outdict['info'] = {}
    outdict['png'] = {}
    try:
        outdict['sta'].update(get_fqc_sta(html_file))
    except Exception as e:
        logging.error(e)
    try:
        with open(html_file,'rt') as h:
            soup = BeautifulSoup(h, 'html.parser')
        try:
            outdict['info'].update(get_fqc_info(soup))
        except Exception as e:
            logging.error(e)
        try:
            outdict['png'].update(get_fqc_png(soup))
        except Exception as e:
            logging.error(e)
            
    except Exception as e:
        logging.error(e)
    return outdict


# %%
def get_fqc_png(soup):
    outdict = {}
    right_info = list(list(soup.body.children)[2])[1:]
    for info in right_info:
        info_tmp = list(info.children)
        key = list(info_tmp[0])[1]
        key = re.sub('\s+','_', key)
        try:
            src = list(info_tmp[1].children)[0]['src']
            src = re.sub('^data:image/png;base64,', '', src)
            value = src
    #         with open(f'{key}.png', 'wb') as h:
    #             imgdata = base64.b64decode(src)
    #             h.write(imgdata)
    #             print(src)
        except Exception as e:
            value = ''
        outdict[key] = value
    return outdict


# %%
def get_fqc_info(soup):
    out_dict = {}
    left_info = list(list(soup.body.children)[1])[1]
    for info in list(left_info.children):
    #     print(info.contents)
        img, title = info.contents
        key = re.sub('\s+', '_', title.string)
        value = re.sub('\]', '', re.sub('\[','', img['alt']))
        out_dict[key] = value
    return out_dict



# %%
def get_fqc_sta(html_file):
    df_sta = pd.read_html(html_file)[0].set_index('Measure')['Value']
    df_sta.index = df_sta.index.str.replace('\s+','_', regex=True).str.replace('%', '')
    return df_sta.to_dict()


# %%
if __name__ == '__main__':
    parse = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parse.add_argument('-i', '--input', required=True, help='fastqc output html file for format')
    parse.add_argument('-o', '--output', default=None, help='output json, default stdout')
    args = parse.parse_args()
    
    logging.basicConfig(level='INFO')
    html_file = args.input
    outdict = get_fqc(html_file)
    if args.output:
        try:
            with open(args.output,'w') as h:
                json.dump(outdict, h, indent=2)
        except Exception as e:
            logging.error(e)
    else:
        print(outdict)

# %%
