import json
import numpy as np
import matplotlib.pyplot as plt
import datetime
from matplotlib import rc, rcParams


dates_query = """{
  "query": {
    "type": "group",
    "nodes": [
      {
        "type": "terminal",
        "label": "text",
        "service": "text",
        "parameters": {
          "operator": "contains_phrase",
          "negation": false,
          "value": "RIBOSOME",
          "attribute": "struct_keywords.pdbx_keywords"
        }
      },
      {
        "type": "terminal",
        "service": "text_chem",
        "parameters": {
          "operator": "in",
          "negation": false,
          "value": [
            "ERY"
          ],
          "attribute": "rcsb_chem_comp_container_identifiers.comp_id"
        }
      }
    ],
    "logical_operator": "and"
  },
  "return_type": "entry",
  "request_options": {
    "scoring_strategy": "combined",
    "sort": [
      {
        "sort_by": "score",
        "direction": "desc"
      }
    ],
    "paginate": {
      "start": 0,
      "rows": 25
    },
    "results_content_type": [
      "experimental"
    ]
  }
}"""
# open the structs.json file and read the data 



with open('custom_report.json', 'r') as infile:
    data = json.load(infile)

dates   = []
heights = []
_       = 0

by_year  = {
    '2000':{
        'xray'        : 0,
        'em'          : 0,
        'emdb_maps'   : 0,
        'xray_resos'  : [],
        'cryoem_resos': []
    }
}

XRAY = 'X-RAY DIFFRACTION'
NMR  = 'SOLUTION NMR'
EM   = 'ELECTRON MICROSCOPY'


for s in data:
    datestr          = str(s['data']['rcsb_accession_info']['initial_release_date']).split("T")[0]
    method           = str(s['data']['exptl'][0]['method']);
    reso             = s['data']['rcsb_entry_info']['resolution_combined']
    maps             = s['data']['rcsb_entry_container_identifiers']['emdb_ids']
    [year,month,day] = datestr.split("-")

    if reso == None :
        continue


    if int(year) < 2000:
        by_year['2000']['emdb_maps'] +=  len (maps) if maps != None else 0

        if method == XRAY:
            by_year['2000']['xray'] += 1
            by_year['2000']['xray_resos'] += reso
        elif method == NMR:
            ...
        elif method == EM:
            by_year['2000']['em'] += 1
            by_year['2000']['cryoem_resos'] += reso
        
    else:

      if year not in by_year:
          by_year[year] = {
          'xray'        : 0,
          'em'          : 0,
          'emdb_maps'   : 0,
          'xray_resos'  : [],
          'cryoem_resos': []
          }

      by_year[year]['emdb_maps'] +=  len (maps) if maps != None else 0
      if method == XRAY:
              by_year[year]['xray'] += 1
              by_year[year]['xray_resos'] += reso

      elif method == NMR:
          ...
      elif method == EM:
              by_year[year]['em'] += 1
              by_year[year]['cryoem_resos'] += reso
    d = datetime.datetime(int(year),int(month),int(day))

    dates.append(d)
    heights.append(_)
    
by_year = sorted(by_year.items())
# ------------------------------------------------------------

d_years       = []
d_ems         = []
d_xrays       = []

d_num_maps    = []
d_cryoem_reso = []
d_xray_reso   = []

for e in by_year:
    _em   = e[1]['em']
    _xray = e[1]['xray']
    _xray_reso = np.mean(e[1]['xray_resos'])
    _cryoem_reso = np.mean(e[1]['cryoem_resos'])
    _num_maps = e[1]['emdb_maps']
    d_years.append(int(e[0]))
    d_ems.append(_em)
    d_xrays.append(_xray)
    d_num_maps.append(_num_maps)
    d_xray_reso.append(_xray_reso)
    d_cryoem_reso.append(_cryoem_reso)

width         = 0.65       # the width of the bars: can also be len(x) sequence
fig, ax = plt.subplots()
rc('axes', linewidth=2)
xraybars = ax.bar(d_years, d_xrays, width,                label='X-RAY DIFFRACTION' , color="black"   , edgecolor="black")
embars   = ax.bar(d_years, d_ems  , width, bottom=d_xrays,label='ELECTRON MICROSCOPY'  ,fill=None , edgecolor="black")
ax.scatter(d_years, d_num_maps, color="blue", marker='*', label="Associated EMDB maps", s=150)
# ax.plot(d_years, d_cryoem_reso, color="cyan", label="Cryo-EM resolution", linewidth=3)
# ax.plot(d_years, d_xray_reso, color="pink", label="X-ray resolution", linewidth=3)
# ax.lines[0].set_linestyle("--")
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
# ax.set_ylabel('Scores')
# ax.set_title('Scores by group and gender')
ax.legend()

print("nume entires", len(d_years))
#  set x axis label to "Year"

ax.set_xlabel('Year'                      , fontsize=24)
# ax.set_ylabel('Number of released RCSB models' , fontsize=24)
# title= "Number of entries in $\it{ribosome.xyz}$ by year"
title= "Number of RIBOSOME structures deposited to $\it{rcsb.org}$"
ax.set_title (title, fontsize=24)
plt.legend(loc=2, prop={'size': 22})


# ax.set_xticklabels(['2000', '2001', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018', '2019'])


for ( bar, i_em, i_xray ) in zip(embars.patches, d_ems, d_xrays):
    value = bar.get_height()
    text = f'{i_em+i_xray}'
    text_x = bar.get_x() + bar.get_width() / 2
    text_y = bar.get_y() + value
    ax.text(text_x, text_y + 3, text, ha='center',color='black',size=18)

fig.set_size_inches(18.5, 10.5)
# plt.show()
fig.savefig('cumulative_entries.png', format='png', dpi=300, bbox_inches='tight')
fig.savefig('cumulative_entries.pdf', format='pdf', dpi=300, bbox_inches='tight')