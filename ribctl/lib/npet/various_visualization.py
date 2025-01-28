from pprint import pprint
from typing import Callable, List, Optional
from matplotlib import pyplot as plt
import pyvista as pv
import json
import numpy as np
import open3d as o3d

from ribctl.lib.npet.kdtree_approach import T

hexcolors = {
    'aliceblue'           : '#F0F8FF',
    'antiquewhite'        : '#FAEBD7',
    'aquamarine'          : '#7FFFD4',
    'azure'               : '#F0FFFF',
    'beige'               : '#F5F5DC',
    'bisque'              : '#FFE4C4',
    'black'               : '#000000',
    'blanchedalmond'      : '#FFEBCD',
    'blue'                : '#0000FF',
    'blueviolet'          : '#8A2BE2',
    'brown'               : '#A52A2A',
    'burlywood'           : '#DEB887',
    'cadetblue'           : '#5F9EA0',
    'chartreuse'          : '#7FFF00',
    'chocolate'           : '#D2691E',
    'coral'               : '#FF7F50',
    'cornflowerblue'      : '#6495ED',
    'cornsilk'            : '#FFF8DC',
    'crimson'             : '#DC143C',
    'cyan'                : '#00FFFF',
    'darkblue'            : '#00008B',
    'darkcyan'            : '#008B8B',
    'darkgoldenrod'       : '#B8860B',
    'darkgray'            : '#A9A9A9',
    'darkgreen'           : '#006400',
    'darkkhaki'           : '#BDB76B',
    'darkmagenta'         : '#8B008B',
    'darkolivegreen'      : '#556B2F',
    'darkorange'          : '#FF8C00',
    'darkorchid'          : '#9932CC',
    'darkred'             : '#8B0000',
    'darksalmon'          : '#E9967A',
    'darkseagreen'        : '#8FBC8F',
    'darkslateblue'       : '#483D8B',
    'darkslategray'       : '#2F4F4F',
    'darkturquoise'       : '#00CED1',
    'darkviolet'          : '#9400D3',
    'deeppink'            : '#FF1493',
    'deepskyblue'         : '#00BFFF',
    'dimgray'             : '#696969',
    'dodgerblue'          : '#1E90FF',
    'firebrick'           : '#B22222',
    'floralwhite'         : '#FFFAF0',
    'forestgreen'         : '#228B22',
    'gainsboro'           : '#DCDCDC',
    'ghostwhite'          : '#F8F8FF',
    'gold'                : '#FFD700',
    'goldenrod'           : '#DAA520',
    'gray'                : '#808080',
    'green'               : '#008000',
    'greenyellow'         : '#ADFF2F',
    'honeydew'            : '#F0FFF0',
    'hotpink'             : '#FF69B4',
    'indianred'           : '#CD5C5C',
    'indigo'              : '#4B0082',
    'ivory'               : '#FFFFF0',
    'khaki'               : '#F0E68C',
    'lavender'            : '#E6E6FA',
    'lavenderblush'       : '#FFF0F5',
    'lawngreen'           : '#7CFC00',
    'lemonchiffon'        : '#FFFACD',
    'lightblue'           : '#ADD8E6',
    'lightcoral'          : '#F08080',
    'lightcyan'           : '#E0FFFF',
    'lightgoldenrodyellow': '#FAFAD2',
    'lightgray'           : '#D3D3D3',
    'lightgreen'          : '#90EE90',
    'lightpink'           : '#FFB6C1',
    'lightsalmon'         : '#FFA07A',
    'lightseagreen'       : '#20B2AA',
    'lightskyblue'        : '#87CEFA',
    'lightslategray'      : '#778899',
    'lightsteelblue'      : '#B0C4DE',
    'lightyellow'         : '#FFFFE0',
    'lime'                : '#00FF00',
    'limegreen'           : '#32CD32',
    'linen'               : '#FAF0E6',
    'magenta'             : '#FF00FF',
    'maroon'              : '#800000',
    'mediumaquamarine'    : '#66CDAA',
    'mediumblue'          : '#0000CD',
    'mediumorchid'        : '#BA55D3',
    'mediumpurple'        : '#9370DB',
    'mediumseagreen'      : '#3CB371',
    'mediumslateblue'     : '#7B68EE',
    'mediumspringgreen'   : '#00FA9A',
    'mediumturquoise'     : '#48D1CC',
    'mediumvioletred'     : '#C71585',
    'midnightblue'        : '#191970',
    'mintcream'           : '#F5FFFA',
    'mistyrose'           : '#FFE4E1',
    'moccasin'            : '#FFE4B5',
    'navajowhite'         : '#FFDEAD',
    'navy'                : '#000080',
    'oldlace'             : '#FDF5E6',
    'olive'               : '#808000',
    'olivedrab'           : '#6B8E23',
    'orange'              : '#FFA500',
    'orangered'           : '#FF4500',
    'orchid'              : '#DA70D6',
    'palegoldenrod'       : '#EEE8AA',
    'palegreen'           : '#98FB98',
    'paleturquoise'       : '#AFEEEE',
    'palevioletred'       : '#DB7093',
    'papayawhip'          : '#FFEFD5',
    'paraview_background' : '#52576e',
    'peachpuff'           : '#FFDAB9',
    'peru'                : '#CD853F',
    'pink'                : '#FFC0CB',
    'plum'                : '#DDA0DD',
    'powderblue'          : '#B0E0E6',
    'purple'              : '#800080',
    'raw_sienna'          : '#965434',
    'rebeccapurple'       : '#663399',
    'red'                 : '#FF0000',
    'rosybrown'           : '#BC8F8F',
    'royalblue'           : '#4169E1',
    'saddlebrown'         : '#8B4513',
    'salmon'              : '#FA8072',
    'sandybrown'          : '#F4A460',
    'seagreen'            : '#2E8B57',
    'seashell'            : '#FFF5EE',
    'sienna'              : '#A0522D',
    'silver'              : '#C0C0C0',
    'skyblue'             : '#87CEEB',
    'slateblue'           : '#6A5ACD',
    'slategray'           : '#708090',
    'snow'                : '#FFFAFA',
    'springgreen'         : '#00FF7F',
    'steelblue'           : '#4682B4',
    'tan'                 : '#D2B48C',
    'teal'                : '#008080',
    'thistle'             : '#D8BFD8',
    'tomato'              : '#FF6347',
    'turquoise'           : '#40E0D0',
    'violet'              : '#EE82EE',
    'wheat'               : '#F5DEB3',
    'white'               : '#FFFFFF',
    'whitesmoke'          : '#F5F5F5',
    'yellow'              : '#FFFF00',
    'yellowgreen'         : '#9ACD32',
    'tab:blue'            : '#1f77b4',
    'tab:orange'          : '#ff7f0e',
    'tab:green'           : '#2ca02c',
    'tab:red'             : '#d62728',
    'tab:purple'          : '#9467bd',
    'tab:brown'           : '#8c564b',
    'tab:pink'            : '#e377c2',
    'tab:gray'            : '#7f7f7f',
    'tab:olive'           : '#bcbd22',
    'tab:cyan'            : '#17becf',
}

diagram_tunnels = {
    "bacteria": [
        "4W29", # GOOD
        "6WD4", # GOOD
        "7UNV",
        "5DM6",

        "5O60",
        "8BUU",
        "6HMA",
        "7RYH", #GOOD

        "7MSZ",
        "7P7T",
        "5MYJ", # GOOD, weird structure
        "7JIL",
        # 3j9y -- good
    ],
    "eukaryota": [
        "6P5N",
        "7QGG",
        "4UG0", # GOOD
        "4U3M",

        "7OYB",
        "7OLC",
        "8EUI",
        "5XXB",

        "4V91",
        "6XU8",
        "4V7E",
        "8P5D",

        "8BTR",
        "3JBO",
        "7CPU",
        "7Q08",
        "6AZ3", #GOOD

        "5T5H",
        "5XY3",
        "7QEP",
    ],
    "archaea": ["4V6U",  # GOOD
                "4V9F",] # GOOD
}

FONT                  = 'courier'
CHAIN_PT_SIZE         = 8
PTC_PT_SIZE           = 40
CHAIN_LANDMARK_COLORS = ["purple","orange", "cornflowerblue", "cornsilk", "crimson", "darkblue", "darkcyan", "darkgoldenrod", "darkgray", "darkgreen", "darkkhaki", "darkmagenta", "darkolivegreen", "darkorange", "darkorchid", "darkred",
"rebeccapurple",
"rosybrown",
"royalblue",
"saddlebrown",
"salmon",
"sandybrown",
"seagreen"]

available_tunnels = {
    "archaea": ["1W2B", "4V6U", "4V9F"],
    "bacteria": [
        "1NKW",
        "1NWX",
        "1NWY",
        "1SM1",
        "1VVJ",
        "1VY4",
        "1VY5",
        "1VY6",
        "1VY7",
        "1XBP",
        "2ZJP",
        "2ZJQ",
        "2ZJR",
        "3CF5",
        "3DLL",
        "3J3V",
        "3J3W",
        "3J7Z",
        "3J9W",
        "3J9Y",
        "3J9Z",
        "3JA1",
        "3JBU",
        "3JCD",
        "3JCJ",
        "3PIO",
        "3PIP",
        "4IO9",
        "4IOA",
        "4IOC",
        "4L47",
        "4L71",
        "4LEL",
        "4LFZ",
        "4LNT",
        "4LSK",
        "4LT8",
        "4P6F",
        "4P70",
        "4TUA",
        "4TUB",
        "4TUC",
        "4TUD",
        "4TUE",
        "4U1U",
        "4U1V",
        "4U20",
        "4U24",
        "4U25",
        "4U26",
        "4U27",
        "4UY8",
        "4V4H",
        "4V4I",
        "4V4J",
        "4V4Q",
        "4V50",
        "4V51",
        "4V52",
        "4V53",
        "4V54",
        "4V56",
        "4V57",
        "4V5A",
        "4V5B",
        "4V5C",
        "4V5D",
        "4V5E",
        "4V5F",
        "4V5G",
        "4V5J",
        "4V5K",
        "4V5L",
        "4V5P",
        "4V5Q",
        "4V5R",
        "4V5S",
        "4V63",
        "4V64",
        "4V67",
        "4V6A",
        "4V6C",
        "4V6D",
        "4V6E",
        "4V6F",
        "4V6G",
        "4V7L",
        "4V7M",
        "4V7P",
        "4V7S",
        "4V7T",
        "4V7U",
        "4V7V",
        "4V7W",
        "4V7X",
        "4V7Y",
        "4V7Z",
        "4V83",
        "4V84",
        "4V85",
        "4V89",
        "4V8A",
        "4V8G",
        "4V8H",
        "4V8I",
        "4V8J",
        "4V8N",
        "4V8O",
        "4V8Q",
        "4V8U",
        "4V8X",
        "4V90",
        "4V95",
        "4V97",
        "4V9A",
        "4V9B",
        "4V9C",
        "4V9D",
        "4V9H",
        "4V9I",
        "4V9J",
        "4V9K",
        "4V9L",
        "4V9N",
        "4V9O",
        "4V9P",
        "4V9Q",
        "4V9R",
        "4V9S",
        "4W29",
        "4W2E",
        "4W2F",
        "4W2G",
        "4W2H",
        "4W2I",
        "4W4G",
        "4WF1",
        "4WOI",
        "4WPO",
        "4WQ1",
        "4WQF",
        "4WQR",
        "4WQU",
        "4WQY",
        "4WR6",
        "4WRA",
        "4WRO",
        "4WSD",
        "4WSM",
        "4WT1",
        "4WT8",
        "4WU1",
        "4WWW",
        "4WZD",
        "4WZO",
        "4XEJ",
        "4Y4O",
        "4Y4P",
        "4YBB",
        "4YPB",
        "4YZV",
        "4Z3S",
        "4Z8C",
        "4ZER",
        "4ZSN",
        "5AFI",
        "5CZP",
        "5DM6",
        "5DM7",
        "5DOX",
        "5DOY",
        "5E7K",
        "5E81",
        "5EL4",
        "5EL5",
        "5EL6",
        "5EL7",
        "5F8K",
        "5FDU",
        "5FDV",
        "5GAD",
        "5GAE",
        "5GAG",
        "5GAH",
        "5H5U",
        "5HAU",
        "5HCP",
        "5HCQ",
        "5HCR",
        "5HD1",
        "5IB8",
        "5IBB",
        "5IMQ",
        "5IT8",
        "5J30",
        "5J3C",
        "5J4B",
        "5J4C",
        "5J4D",
        "5J5B",
        "5J7L",
        "5J88",
        "5J8A",
        "5J8B",
        "5J91",
        "5JC9",
        "5JTE",
        "5JU8",
        "5KCR",
        "5KCS",
        "5KPS",
        "5KPW",
        "5KPX",
        "5L3P",
        "5LZA",
        "5LZD",
        "5LZE",
        "5MDV",
        "5MDW",
        "5MDY",
        "5MDZ",
        "5MGP",
        "5MYJ",
        "5NDJ",
        "5NDK",
        "5NGM",
        "5NP6",
        "5NWY",
        "5O2R",
        "5O60",
        "5OT7",
        "5T7V",
        "5U4I",
        "5U9F",
        "5U9G",
        "5UQ7",
        "5UQ8",
        "5UYK",
        "5UYL",
        "5UYM",
        "5UYP",
        "5UYQ",
        "5V7Q",
        "5V8I",
        "5V93",
        "5VP2",
        "5VPO",
        "5VPP",
        "5W4K",
        "5WDT",
        "5WE4",
        "5WE6",
        "5WF0",
        "5WFK",
        "5WFS",
        "5WIS",
        "5WIT",
        "5XYM",
        "5ZEB",
        "5ZEP",
        "5ZET",
        "5ZLU",
        "6B4V",
        "6BOH",
        "6BOK",
        "6BU8",
        "6BUW",
        "6BY1",
        "6BZ6",
        "6BZ7",
        "6BZ8",
        "6C4I",
        "6C5L",
        "6CAE",
        "6CFJ",
        "6CFK",
        "6CFL",
        "6CZR",
        "6DNC",
        "6DZI",
        "6DZP",
        "6ENF",
        "6ENJ",
        "6ENU",
        "6FKR",
        "6GC8",
        "6GSJ",
        "6GSK",
        "6GSL",
        "6GWT",
        "6GXM",
        "6GXN",
        "6GXO",
        "6GZQ",
        "6HA1",
        "6HA8",
        "6HMA",
        "6HTQ",
        "6I0Y",
        "6I7V",
        "6N1D",
        "6N9E",
        "6N9F",
        "6ND5",
        "6ND6",
        "6NDK",
        "6NSH",
        "6NTA",
        "6NUO",
        "6NWY",
        "6O3M",
        "6O97",
        "6O9J",
        "6OF1",
        "6OF6",
        "6OFX",
        "6OG7",
        "6OGF",
        "6OGI",
        "6OJ2",
        "6OM6",
        "6OPE",
        "6ORD",
        "6ORE",
        "6ORL",
        "6OSK",
        "6OSQ",
        "6OT3",
        "6OTR",
        "6OUO",
        "6OXA",
        "6OXI",
        "6PJ6",
        "6Q95",
        "6Q97",
        "6Q9A",
        "6QNQ",
        "6QNR",
        "6S0K",
        "6S0X",
        "6S0Z",
        "6S12",
        "6S13",
        "6SZS",
        "6TBV",
        "6TC3",
        "6U48",
        "6UCQ",
        "6UO1",
        "6V39",
        "6V3A",
        "6V3B",
        "6V3D",
        "6VU3",
        "6VWL",
        "6VWM",
        "6VWN",
        "6VYQ",
        "6VYR",
        "6VYS",
        "6WD0",
        "6WD1",
        "6WD2",
        "6WD3",
        "6WD4",
        "6WD5",
        "6WD6",
        "6WD7",
        "6WD8",
        "6WD9",
        "6WDA",
        "6WDD",
        "6WDE",
        "6WDF",
        "6WDG",
        "6WDJ",
        "6WDK",
        "6WDL",
        "6WDM",
        "6WRU",
        "6X7F",
        "6X7K",
        "6XDQ",
        "6XHV",
        "6XHW",
        "6XHX",
        "6XHY",
        "6XQD",
        "6XQE",
        "6XZ7",
        "6XZA",
        "6XZB",
        "6Y69",
        "6YSR",
        "6YSS",
        "6YST",
        "6YSU",
        "7AQC",
        "7AQD",
        "7ASO",
        "7ASP",
        "7AZO",
        "7AZS",
        "7B5K",
        "7BL2",
        "7BL3",
        "7BL4",
        "7BV8",
        "7D6Z",
        "7F0D",
        "7JQL",
        "7JQM",
        "7JSS",
        "7JSW",
        "7JSZ",
        "7JT1",
        "7JT2",
        "7JT3",
        "7K00",
        "7K50",
        "7K51",
        "7K52",
        "7K53",
        "7K54",
        "7K55",
        "7KGB",
        "7LH5",
        "7LV0",
        "7LVK",
        "7M4V",
        "7M4W",
        "7M4X",
        "7M4Y",
        "7M4Z",
        "7MD7",
        "7MSC",
        "7MSH",
        "7MSM",
        "7MSZ",
        "7MT2",
        "7MT3",
        "7MT7",
        "7N1P",
        "7N2C",
        "7N2U",
        "7N2V",
        "7N30",
        "7N31",
        "7NSO",
        "7NSP",
        "7NSQ",
        "7NWT",
        "7NWW",
        "7O5B",
        "7OIF",
        "7OIG",
        "7OII",
        "7OIZ",
        "7OJ0",
        "7OPE",
        "7OT5",
        "7P3K",
        "7P7T",
        "7PJS",
        "7PJV",
        "7PJY",
        "7QG8",
        "7QGN",
        "7QGU",
        "7QH4",
        "7QV1",
        "7QV3",
        "7RQ8",
        "7RQ9",
        "7RQA",
        "7RQB",
        "7RQC",
        "7RQD",
        "7RQE",
        "7RYF",
        "7RYG",
        "7RYH",
        "7S0S",
        "7S1G",
        "7S1H",
        "7S1I",
        "7S1J",
        "7S1K",
        "7SFR",
        "7SS9",
        "7SSD",
        "7SSL",
        "7SSN",
        "7SSO",
        "7SSW",
        "7ST2",
        "7ST6",
        "7ST7",
        "7TOS",
        "7U2H",
        "7U2I",
        "7U2J",
        "7UG7",
        "7UNR",
        "7UNU",
        "7UNV",
        "7UNW",
        "7XAM",
        "7Y41",
        "8BUU",
    ],
    "eukaryota": [
        "3J6X",
        "3J6Y",
        "3J77",
        "3J78",
        "3J79",
        "3J7O",
        "3J7P",
        "3J7Q",
        "3J7R",
        "3J92",
        "3JAG",
        "3JAH",
        "3JAI",
        "3JAJ",
        "3JAN",
        "3JBN",
        "3JBO",
        "3JBP",
        "3JCS",
        "3JCT",
        "4D5Y",
        "4D67",
        "4U3M",
        "4U3N",
        "4U3U",
        "4U4N",
        "4U4O",
        "4U4Q",
        "4U4R",
        "4U4U",
        "4U4Y",
        "4U4Z",
        "4U50",
        "4U51",
        "4U52",
        "4U53",
        "4U55",
        "4U56",
        "4U6F",
        "4UG0",
        "4UJC",
        "4UJD",
        "4UJE",
        "4V3P",
        "4V6I",
        "4V6W",
        "4V6X",
        "4V7E",
        "4V7H",
        "4V7R",
        "4V88",
        "4V8M",
        "4V8P",
        "4V8T",
        "4V8Y",
        "4V8Z",
        "4V91",
        "5AJ0",
        "5APN",
        "5APO",
        "5DAT",
        "5DC3",
        "5DGE",
        "5DGF",
        "5DGV",
        "5FCI",
        "5FCJ",
        "5GAK",
        "5H4P",
        "5I4L",
        "5IT7",
        "5JCS",
        "5JUO",
        "5JUP",
        "5JUS",
        "5JUT",
        "5JUU",
        "5LKS",
        "5LYB",
        "5LZS",
        "5LZT",
        "5LZU",
        "5LZV",
        "5LZW",
        "5LZX",
        "5LZY",
        "5LZZ",
        "5M1J",
        "5MC6",
        "5MEI",
        "5NDG",
        "5NDV",
        "5NDW",
        "5OBM",
        "5ON6",
        "5T2A",
        "5T2C",
        "5T5H",
        "5T62",
        "5T6R",
        "5TBW",
        "5TGA",
        "5TGM",
        "5UMD",
        "5XXB",
        "5XY3",
        "6AZ3",
        "6D90",
        "6D9J",
        "6GQ1",
        "6GQB",
        "6GQV",
        "6GZ3",
        "6GZ4",
        "6GZ5",
        "6HCF",
        "6HCJ",
        "6HCM",
        "6HCQ",
        "6HHQ",
        "6IP5",
        "6IP6",
        "6IP8",
        "6M62",
        "6N8J",
        "6N8K",
        "6N8L",
        "6N8M",
        "6N8N",
        "6N8O",
        "6OIG",
        "6OLE",
        "6OLF",
        "6OLG",
        "6OLI",
        "6OLZ",
        "6OM0",
        "6OM7",
        "6P5I",
        "6P5J",
        "6P5K",
        "6P5N",
        "6QIK",
        "6QT0",
        "6QTZ",
        "6R5Q",
        "6R6G",
        "6R6P",
        "6R7Q",
        "6R84",
        "6R86",
        "6R87",
        "6RI5",
        "6RM3",
        "6RZZ",
        "6S05",
        "6SGC",
        "6SWA",
        "6T59",
        "6UZ7",
        "6W6L",
        "6WOO",
        "6XIQ",
        "6XIR",
        "6XU6",
        "6XU7",
        "6XU8",
        "6Y0G",
        "6Y2L",
        "6Y57",
        "6YLX",
        "6Z6J",
        "6Z6K",
        "6Z6L",
        "6Z6M",
        "6Z6N",
        "6ZU5",
        "7A01",
        "7AZY",
        "7BHP",
        "7CPU",
        "7CPV",
        "7LS1",
        "7LS2",
        "7MDZ",
        "7NFX",
        "7NRC",
        "7NRD",
        "7NWG",
        "7NWH",
        "7NWI",
        "7OF1",
        "7OH3",
        "7OHQ",
        "7OLC",
        "7OLD",
        "7OSA",
        "7OYA",
        "7OYB",
        "7OYC",
        "7PWG",
        "7PWO",
        "7PZY",
        "7Q08",
        "7Q0F",
        "7Q0P",
        "7Q0R",
        "7QCA",
        "7QEP",
        "7QGG",
        "7QWQ",
        "7QWR",
        "7QWS",
        "7RR5",
        "7TM3",
        "7TOP",
        "7TOQ",
        "7TOR",
        "7TUT",
        "7UG6",
        "7Z34",
        "7ZPQ",
        "7ZRS",
        "7ZS5",
        "7ZUW",
        "7ZUX",
        "7ZW0",
        "8BIP",
        "8BJQ",
        "8BQD",
        "8BQX",
        "8BR8",
        "8BRM",
        "8BSI",
        "8BSJ",
        "8BTD",
        "8BTR",
        "8EUG",
        "8EUI",
        "8HFR",
        "8IR3",
        "8P5D",
        "8P60",
        "8Q5I",
    ],
}

dbscan_pairs = [
    (2,33),
    (2.3, 57),
    (3, 123),
    (3.2, 145),
    (3.5, 175),
    (4.2, 280),
    (5, 490),
    (5.5,600)
]


def visualize_DBSCAN_CLUSTERS_particular_eps_minnbrs(dbscan_cluster_dict: dict[int, list], eps, min_nbrs, base_point, axis_point, refined, radius, height):
    plotter = pv.Plotter()
    plotter.subplot(0, 0)

    # Print cluster sizes
    cluster_sizes = {k: len(v) for k, v in dbscan_cluster_dict.items()}
    for k, v in cluster_sizes.items():
        print(f"Cluster {k} has {v} points.")

    # Find the largest non-noise cluster
    largest_cluster = max(
        ((k, len(v)) for k, v in dbscan_cluster_dict.items() if k != -1),
        key=lambda x: x[1]
    )[0]

    # Generate color palette for regular clusters
    clusters_palette = dict(zip(range(-1, 60), plt.cm.terrain(np.linspace(0, 1, 60))))
    for k, v in clusters_palette.items():
        clusters_palette[k] = [*v[:3], 0.5]

    # Create cylinder
    # Calculate direction vector from base_point to axis_point
    direction = np.array(axis_point) - np.array(base_point)
    direction = direction / np.linalg.norm(direction)
    center = np.array(base_point) + (direction * height/2)
    
    # Create cylinder jcentered at base_point
    cylinder = pv.Cylinder(
        center=center,  # Use calculated center
        direction=direction,  # Remove the negative sign
        radius=radius,
        height=height
    )
    
    # Add cylinder to plot with transparency
    plotter.add_mesh(cylinder, opacity=0.1, color='gray', label='Cylinder', style='wireframe')

    # Add each cluster separately
    for cluster_label, coordinates in dbscan_cluster_dict.items():
        points = np.array(coordinates)
        if cluster_label == -1:
            # Noise points
            # ...
            plotter.add_points(
                points,
                color      = 'gray',
                opacity    = 0.1,
                point_size = 1,
                label      = 'Noise',
                style      = 'points_gaussian',
                emissive   = True
            )

        elif cluster_label == largest_cluster:
            # Largest cluster - blue spheres
            print(points.shape)
            plotter.add_points(
                points,
                color='cyan',
                point_size=5,
                opacity=0.15,
                label=f'Cluster {cluster_label} (Largest)',
                style='points',
                render_points_as_spheres=True
            )
        else:
            # Regular clusters
            color = clusters_palette[(cluster_label * 2) % len(clusters_palette)]
            plotter.add_points(
                points,
                # color=color[:3],
                color='gray',
                opacity=0.3,
                point_size=2,
                label=f'Cluster {cluster_label}',
                style='points'
            )

    plotter.add_points(
        refined,
        color='blue',
        point_size=3,
        opacity=1,
        label='refined_points',
        style='points',
        render_points_as_spheres=True
    )
    # Add reference points
    plotter.add_points(
        np.array([base_point]),
        color='red',
        point_size=20,
        label='Base Point',
        style='points',
        render_points_as_spheres=True
    )
    plotter.add_points(
        np.array([axis_point]),
        color='red',
        point_size=20,
        label='Axis Point',
        style='points',
        render_points_as_spheres=True
    )

    # Add text information
    plotter.add_text(
        f'eps: {eps}\nmin_nbrs: {min_nbrs}',
        position='upper_left',
        font_size=20,
        shadow=True,
        font=FONT,
        color='black'
    )

    plotter.show()

def DBSCAN_CLUSTERS_visualize_largest(positive_space: np.ndarray, dbscan_cluster_dict: dict[int, list], selected_cluster: np.ndarray, gif:bool=False, gif_name:str|None=None):
    plotter               = pv.Plotter(shape=(1, 2), off_screen=True)
    plotter.subplot(0,0)
    n_labels = 7
    plotter.add_axes(line_width=2,cone_radius=0.3, shaft_length=2, tip_length=1, ambient=1, label_size=(0.2, 0.6))
    plotter.show_grid( n_xlabels=n_labels, n_ylabels=n_labels, n_zlabels=n_labels, font_size = 8)
    #? Visualize all clusters
    for k, v in dbscan_cluster_dict.items():
        print("Cluster {} has {} points.".format(k, len(v)))
    clusters_palette = dict(zip(range(-1, 60), plt.cm.terrain(np.linspace(0, 1, 60))))
    for k, v in clusters_palette.items():
        clusters_palette[k] = [*v[:3], 0.5]
    combined_cluster_colors = []
    combined_cluster_points = []

    for dbscan_label, coordinates in dbscan_cluster_dict.items():
        combined_cluster_points.extend(coordinates)
        combined_cluster_colors.extend([clusters_palette[( dbscan_label * 5 )%len(clusters_palette)]   if dbscan_label != -1 else [0, 0, 0, 0.1]] * len(coordinates) )

    ptcloud_all_clusters         = pv.PolyData(combined_cluster_points)
    ptcloud_all_clusters["rgba"] = combined_cluster_colors

    plotter.add_mesh(ptcloud_all_clusters, scalars="rgba", rgb=True, show_scalar_bar=False)

    # ? Visualize selected cluster
    plotter.subplot(0,1)

    lc = selected_cluster
    print("Max vals in selected cluster:", [[np.min(lc[:,0]), np.max(lc[:,0])], [np.min(lc[:,1]), np.max(lc[:,1])],[np.min(lc[:,2]), np.max(lc[:,2])] ])

    n_labels = 7
    plotter.add_axes(line_width=2,cone_radius=0.3, shaft_length=2, tip_length=1, ambient=1, label_size=(0.2, 0.6))
    plotter.show_grid( n_xlabels=n_labels, n_ylabels=n_labels, n_zlabels=n_labels, font_size = 12)

    rgbas_cluster       = [[15, 10, 221, 1] for datapoint in selected_cluster]
    rgbas_positive      = np.array([[205, 209, 228, 0.2] for _ in positive_space])
    combined            = np.concatenate([selected_cluster, positive_space])
    rgbas_combined      = np.concatenate([rgbas_cluster, rgbas_positive])

    # list(selected_cluster).extend(list(positive_space))
    # list(rgbas_cluster).extend(list(rgbas_positive))
    # combined         = selected_cluster.extend(positive_space)
    # rgbas_combined   = np.concatenate([rgbas_cluster, rgbas_positie])
    container_points = [*list(selected_cluster),*list(positive_space)]
    # container_rgbas  = [*list(rgbas_cluster),*list(rgbas_positive)]

    np.save("selected_cluster.npy", selected_cluster)
    np.save("positive_space.npy", positive_space)
    point_cloud         = pv.PolyData(selected_cluster)
    point_cloud["rgba"] = rgbas_cluster

    point_cloud         = pv.PolyData(container_points)
    # point_cloud["rgba"] = container_rgbas
    plotter.add_points(point_cloud, rgb=True, show_scalar_bar=True)
    if gif:
        output_gif = gif_name
        # plotter.camera.zoom(1.5)
        plotter.open_gif(output_gif)

        # Rotate the camera 360 degrees
        for angle in range(0, 360, 2):  # 5 degree steps
            plotter.camera.azimuth = angle
            plotter.write_frame()
        plotter.close()
        print(f"GIF saved as {output_gif}")
    else:
        plotter.show()

def visualize_mesh(mesh_path, rcsb_id:str|None=None, gif:bool=False, gif_name:str|None=None):
    plotter = pv.Plotter(off_screen=gif)


    _ = plotter.add_mesh(pv.read(mesh_path), opacity=0.8)
    plotter.add_axes(line_width=2,cone_radius=0.7, shaft_length=0.7, tip_length=0.3, ambient=0.5, label_size=(0.2, 0.8))
    plotter.add_text('RCSB_ID:{}'.format(rcsb_id if rcsb_id is not None else "" ), position='upper_right', font_size=14, shadow=True, font='courier', color='black')
    plotter.show_grid( n_xlabels=8, n_ylabels=8, n_zlabels=8, font_size = 8)

    if gif:
        output_gif = gif_name
        # plotter.camera.zoom(1.5)
        plotter.open_gif(output_gif)

        # Rotate the camera 360 degrees
        for angle in range(0, 360, 2):  # 5 degree steps
            plotter.camera.azimuth = angle
            plotter.write_frame()
        plotter.close()
        print(f"GIF saved as {output_gif}")
    else:
        plotter.show()

def visualize_pointcloud(ptcloud,  rcsb_id:str|None=None, gif:bool=False, gif_name:str|None=None):
    plotter              = pv.Plotter(off_screen=gif)
    n_labels = 7
    plotter.show_grid( n_xlabels=n_labels, n_ylabels=n_labels, n_zlabels=n_labels, font_size = 8)
    plotter.add_axes(line_width=2,cone_radius=0.3, shaft_length=2, tip_length=1, ambient=1, label_size=(0.2, 0.6))
    plotter.add_text('RCSB_ID:{}'.format(rcsb_id if rcsb_id is not None else "" ), position='upper_right', font_size=14, shadow=True, font='courier', color='black')
    plotter.add_points(pv.PolyData(ptcloud), color='black', point_size=3,  opacity=0.3)
    plotter.show()

def visualize_pointcloud_axis(ptcloud,aux_ptcloud, base_point, axis_point, radius=0.1, height=None, rcsb_id:str|None=None, gif:bool=False, gif_name:str|None=None):
    # Convert inputs to numpy arrays to ensure compatibility
    base_point = np.asarray(base_point)
    axis_point = np.asarray(axis_point)
    ptcloud = np.asarray(ptcloud)
    
    # Calculate the original axis vector
    original_axis = axis_point - base_point
    original_axis_length = np.linalg.norm(original_axis)
    original_axis = original_axis / original_axis_length
    
    # Set height to axis length if not specified
    if height is None:
        height = original_axis_length
    
    # Define the target z-axis (pointing down)
    target_z_axis = np.array([0, 0, -1])
    
    # Calculate the rotation matrix to align the original axis with the z-axis
    rotation_axis = np.cross(original_axis, target_z_axis)
    
    # Calculate the angle between the vectors
    angle = np.arccos(np.dot(original_axis, target_z_axis))
    
    # Create rotation matrix using Rodrigues' rotation formula
    if np.linalg.norm(rotation_axis) > 1e-10:
        rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
        K = np.array([
            [0, -rotation_axis[2], rotation_axis[1]],
            [rotation_axis[2], 0, -rotation_axis[0]],
            [-rotation_axis[1], rotation_axis[0], 0]
        ])
        R = (np.eye(3) + np.sin(angle) * K + 
             (1 - np.cos(angle)) * np.dot(K, K))
    else:
        # If the vectors are already aligned, use identity matrix
        R = np.eye(3)
    
    # Transform the point cloud
    transformed_ptcloud = np.dot(ptcloud - base_point, R.T)
    transformed_aux_pts = np.dot(aux_ptcloud - base_point, R.T)
    
    # Transform the base and axis points
    transformed_base_point = np.dot(base_point - base_point, R.T)
    transformed_axis_point = np.dot(axis_point - base_point, R.T)
    
    # Create the plotter
    plotter = pv.Plotter()
    
    # Add visualization elements
    n_labels = 7
    # plotter.show_grid(n_xlabels=n_labels, n_ylabels=n_labels, n_zlabels=n_labels, font_size=8)
    # plotter.add_axes(line_width=2, cone_radius=0.3, shaft_length=2, tip_length=1, ambient=1, label_size=(0.2, 0.6))
    
    # # Add RCSB ID text
    # plotter.add_text(f'RCSB_ID:{rcsb_id if rcsb_id is not None else ""}', 
    #                  position='upper_right', font_size=14, 
    #                  shadow=True, font='courier', color='black')
    
    # Create and add transformed cylinder
    cylinder = pv.Cylinder(
        center=transformed_base_point + (height/2) * np.array([0, 0, -1]),
        direction=[0, 0, -1],
        radius=radius,
        height=height,
        resolution=50
    ).triangulate()
    
    # Add cylinder with visible edges
    # plotter.add_mesh(cylinder, color='blue', opacity=0.3, show_edges=True, edge_color='black')

    plotter.add_mesh( cylinder, style='surface', line_width=3, show_edges=True, opacity=0.05, color='lightgreen', silhouette=True )
    plotter.add_points(transformed_base_point, color='red', point_size=12, render_points_as_spheres=True)
    plotter.add_points(transformed_axis_point, color='red', point_size=12, render_points_as_spheres=True)
    plotter.add_points(pv.PolyData(transformed_ptcloud), render_points_as_spheres=True,  style='points', point_size=3, opacity=0.4)
    plotter.add_points(pv.PolyData(transformed_aux_pts), render_points_as_spheres=True,  style='points', color='blue', point_size=5, opacity=0.8)
    
    # Show the plotr
    plotter.show()
    
    camera_position = plotter.camera_position
    # camera_focal_point = plotter.camera_focal_point
    # camera_viewup = plotter.camera_viewup

    # Then in a separate render
    plotter_render = pv.Plotter(off_screen=True)
    plotter_render.add_mesh( cylinder, style='surface', line_width=3, show_edges=True, opacity=0.05, color='lightgreen', silhouette=True )
    plotter_render.add_points(transformed_base_point, color='red', point_size=12, render_points_as_spheres=True)
    plotter_render.add_points(transformed_axis_point, color='red', point_size=12, render_points_as_spheres=True)
    plotter_render.add_points(pv.PolyData(transformed_ptcloud), render_points_as_spheres=True,  style='points', point_size=3, opacity=0.4)
    plotter_render.add_points(pv.PolyData(transformed_aux_pts), render_points_as_spheres=True,  style='points', color='blue', point_size=5, opacity=0.8)
    
    plotter_render.camera_position = camera_position
    plotter_render.screenshot('4ug0_ptcloud.png', transparent_background=True, window_size=(1920,1280))
    
    return transformed_ptcloud

def create_cylinder_from_points(base_point, axis_point, radius, height):
    """
    Create a cylinder using two points: a base point and a point defining the axis direction.
    
    Parameters:
    base_point: array-like, [x, y, z] starting point of cylinder
    axis_point: array-like, [x, y, z] point defining cylinder axis direction
    radius: float, radius of cylinder
    height: float, total height of cylinder (can extend past axis_point)
    
    Returns:
    cylinder: PyVista cylinder mesh (triangulated)
    direction: numpy array of normalized direction vector
    """
    base_point = np.array(base_point)
    axis_point = np.array(axis_point)
    
    direction = axis_point - base_point
    direction = direction / np.linalg.norm(direction)
    
    cylinder = pv.Cylinder(
        center=base_point + (height/2) * direction,
        direction=direction,
        radius=radius,
        height=height,
        resolution=30
    ).triangulate()
    
    return cylinder, direction

def visualize_pointclouds(ptcloud1:np.ndarray, ptcloud2:np.ndarray, background_positive:np.ndarray, gif:bool=False, gif_name:str|None=None):
    plotter               = pv.Plotter(shape=(1, 2), off_screen=gif)
    plotter.subplot(0,0)
    n_labels = 7
    plotter.add_axes(line_width=2,cone_radius=0.3, shaft_length=2, tip_length=1, ambient=1, label_size=(0.2, 0.6))
    plotter.show_grid( n_xlabels=n_labels, n_ylabels=n_labels, n_zlabels=n_labels, font_size = 12)

    rgbas_cluster1       = [[15, 100, 21, 1] for datapoint in ptcloud1]
    rgbas_positive1      = np.array([[205, 209, 228, 0.2] for _ in background_positive])

    combined1            = np.concatenate([ptcloud1, background_positive])
    rgbas_combined1      = np.concatenate([rgbas_cluster1, rgbas_positive1])

    point_cloud1         = pv.PolyData(combined1)
    point_cloud1["rgba"] = rgbas_combined1

    plotter.add_points(point_cloud1, scalars="rgba", rgb=True, show_scalar_bar=False)


    # ? Visualize selected cluster
    plotter.subplot(0,1)
    n_labels = 7
    plotter.add_axes(line_width=2,cone_radius=0.3, shaft_length=2, tip_length=1, ambient=1, label_size=(0.2, 0.6))
    plotter.show_grid( n_xlabels=n_labels, n_ylabels=n_labels, n_zlabels=n_labels, font_size = 12)

    rgbas_cluster2       = [[15, 10, 221, 1] for datapoint in ptcloud2]
    rgbas_positive2      = np.array([[205, 160, 200, 0.2] for _ in background_positive])
    combined2            = np.concatenate([ptcloud2, background_positive])
    rgbas_combined2      = np.concatenate([rgbas_cluster2, rgbas_positive2])
    point_cloud2         = pv.PolyData(combined2)
    point_cloud2["rgba"] = rgbas_combined2
    plotter.add_points(point_cloud2, scalars="rgba", rgb=True, show_scalar_bar=False)


    if gif:
        output_gif = gif_name
        # plotter.camera.zoom(1.5)
        plotter.open_gif(output_gif)

        # Rotate the camera 360 degrees
        for angle in range(0, 360, 2):  # 5 degree steps
            plotter.camera.azimuth = angle
            plotter.write_frame()
        plotter.close()
        print(f"GIF saved as {output_gif}")
    else:
        plotter.show()


    plotter.show()

def visualize_clipping_result(
    original_points: np.ndarray,
    clipped_points: np.ndarray,
    mesh_path: str,
    show_mesh: bool = True,
):
    """
    Visualizes the original points, clipped points, and the clipping mesh.

    Parameters:
        original_points (np.ndarray): The original point cloud
        clipped_points (np.ndarray): The clipped point cloud
        mesh_path (str): Path to the mesh file
        show_mesh (bool): Whether to show the mesh or not
    """
    p = pv.Plotter()

    # Add original points in red
    original_cloud = pv.PolyData(original_points)
    p.add_mesh(
        original_cloud,
        color="red",
        point_size=5,
        render_points_as_spheres=True,
        label="Original Points",
    )

    # Add clipped points in blue
    clipped_cloud = pv.PolyData(clipped_points)
    p.add_mesh(
        clipped_cloud,
        color="blue",
        point_size=5,
        render_points_as_spheres=True,
        label="Clipped Points",
    )

    # Add mesh if requested
    if show_mesh:
        mesh = pv.read(mesh_path)
        p.add_mesh(
            mesh, style="wireframe", color="gray", opacity=0.5, label="Clipping Mesh"
        )

    p.add_legend()
    p.show()

def visualize_filtered_residues(
    filtered_residues: List[T],
    all_residues: Optional[List[T]],  # If None, won't show unfiltered residues
    base_point: np.ndarray,
    axis_point: np.ndarray,
    radius: float,
    height: float,
    position_getter: Callable[[T], np.ndarray] = lambda x: x.center_of_mass(),
    point_size: float = 5,
    opacity: float = 0.3,
    show: bool = True,
    screenshot_path: Optional[str] = None,
    window_size: tuple = (1024, 768)
) -> Optional[pv.Plotter]:

    """
    Visualize filtered residues alongside the cylinder that was used for filtering.
    
    Parameters:
    -----------
    filtered_residues : List[T]
        List of residues that passed the cylinder filter
    all_residues : Optional[List[T]]
        Complete list of residues before filtering. If provided, will show
        unfiltered residues in gray
    base_point : np.ndarray
        Center point of cylinder base
    axis_point : np.ndarray
        Point defining cylinder axis direction
    radius : float
        Radius of cylinder
    height : float
        Height of cylinder
    position_getter : Callable[[T], np.ndarray]
        Function to extract position from residue object
    point_size : float
        Size of points in visualization
    opacity : float
        Opacity of cylinder (0.0-1.0)
    show : bool
        Whether to display the plot immediately
    screenshot_path : Optional[str]
        If provided, saves screenshot to this path
    window_size : tuple
        Size of the visualization window (width, height)
        
    Returns:
    --------
    Optional[pv.Plotter]
        Returns plotter if show=False, None otherwise
    """
    # Initialize plotter
    plotter = pv.Plotter(window_size=window_size)
    
    # Get positions of filtered residues
    filtered_positions = np.array([
        position_getter(res) for res in filtered_residues
    ])
    
    # Create and add cylinder
    direction = axis_point - base_point
    direction = direction / np.linalg.norm(direction)
    
    cylinder = pv.Cylinder(
        center=base_point + (height/2) * direction,
        direction=direction,
        radius=radius,
        height=height,
        resolution=30
    )
    
    # Add cylinder with wireframe and solid style
    plotter.add_mesh(
        cylinder, 
        style='wireframe', 
        color='red', 
        line_width=2, 
        label='Cylinder'
    )
    plotter.add_mesh(
        cylinder, 
        style='surface', 
        color='red', 
        opacity=opacity,
    )
    
    # If all_residues provided, show unfiltered residues in gray
    if all_residues is not None:
        all_positions = np.array([
            position_getter(res) for res in all_residues
        ])
        plotter.add_points(
            all_positions, 
            color='gray', 
            point_size=point_size-2, 
            opacity=0.3,
            label='Unfiltered Residues'
        )
    
    # Add filtered points
    plotter.add_points(
        filtered_positions, 
        color='blue', 
        point_size=point_size,
        label='Filtered Residues'
    )
    
    # Add base and axis points for reference
    plotter.add_points(
        np.array([base_point]), 
        color='green', 
        point_size=point_size*2, 
        label='Base Point'
    )
    plotter.add_points(
        np.array([axis_point]), 
        color='yellow', 
        point_size=point_size*2, 
        label='Axis Point'
    )
    
    # Add axis line
    axis_line = pv.Line(base_point, axis_point)
    plotter.add_mesh(
        axis_line, 
        color='white', 
        line_width=2, 
        label='Cylinder Axis'
    )
    
    # Add legend
    plotter.add_legend()
    
    # Set camera position for better initial view
    plotter.camera_position = 'iso'
    plotter.enable_eye_dome_lighting()  # Improves point visibility
    
    # Add text with stats
    stats_text = (
        f'Total filtered residues: {len(filtered_residues)}\n'
        f'Cylinder height: {height:.1f}\n'
        f'Cylinder radius: {radius:.1f}'
    )
    if all_residues:
        stats_text = f'Total residues: {len(all_residues)}\n' + stats_text
        
    plotter.add_text(
        stats_text,
        position='upper_left',
        font_size=12,
        color='white'
    )
        
    if show:
        plotter.show()
        return None
        
    return plotter

def ptcloud_convex_hull_points_and_gif(
    pointcloud: np.ndarray,
    ALPHA: float,
    TOLERANCE: float,
    output_path: str,
    n_frames: int = 180,
) -> np.ndarray:

    from tqdm import tqdm
    from PIL import Image

    assert pointcloud is not None

    # Create the convex hull
    cloud = pv.PolyData(pointcloud)
    grid = cloud.delaunay_3d(alpha=ALPHA, tol=TOLERANCE, offset=2, progress_bar=True)
    convex_hull = grid.extract_surface()

    # Set up the plotter
    plotter = pv.Plotter(off_screen=True)
    plotter.add_mesh(convex_hull, color="lightblue", show_edges=True)
    plotter.add_points(cloud, color="red", point_size=5)

    # Set up the camera
    plotter.camera_position = "xy"

    # Create frames
    frames = []
    for i in tqdm(range(n_frames)):
        plotter.camera.azimuth = i * (360 / n_frames)
        plotter.render()
        image = plotter.screenshot(transparent_background=False, return_img=True)
        frames.append(Image.fromarray(image))

    # Save as GIF
    frames[0].save(
        output_path, save_all=True, append_images=frames[1:], duration=50, loop=0
    )

    print(f"GIF saved to {output_path}")

    return convex_hull.points

def estimate_normals_and_create_gif(
    convex_hull_surface_pts: np.ndarray,
    output_path: str,
    kdtree_radius=None,
    kdtree_max_nn=None,
    correction_tangent_planes_n=None,
):
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(convex_hull_surface_pts)
    pcd.estimate_normals(
        search_param=o3d.geometry.KDTreeSearchParamHybrid(
            radius=kdtree_radius, max_nn=kdtree_max_nn
        )
    )
    pcd.orient_normals_consistent_tangent_plane(k=correction_tangent_planes_n)

    # Set up the visualizer
    vis = o3d.visualization.Visualizer()
    vis.create_window()
    vis.add_geometry(pcd)
    opt = vis.get_render_option()
    opt.point_size = 5
    opt.point_show_normal = True

    # Set up the camera
    ctr = vis.get_view_control()
    ctr.set_zoom(0.8)

    # Capture frames
    frames = []
    for i in range(360):
        ctr.rotate(10.0, 0.0)  # Rotate 10 degrees around the z-axis
        vis.update_geometry(pcd)
        vis.poll_events()
        vis.update_renderer()
        image = vis.capture_screen_float_buffer(False)
        frames.append(Image.fromarray((np.array(image) * 255).astype(np.uint8)))

    vis.destroy_window()

    frames[0].save(
        output_path, save_all=True, append_images=frames[1:], duration=50, loop=0
    )

    print(f"GIF saved to {output_path}")
