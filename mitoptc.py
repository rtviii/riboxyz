from typing import TypeVar


rcsb_ids = [
      "2FTC",
      "3J6B",
      "3J7Y",
      "3J9M",
      "3JD5",
      "4CE4",
      "4V19",
      "5AJ3",
      "5AJ4",
      "5MRC",
      "5MRE",
      "5MRF",
      "5OOL",
      "5OOM",
      "6GAW",
      "6GAZ",
      "6GB2",
      "6I9R",
      "6NU2",
      "6NU3",
      "6RW4",
      "6RW5",
      "6VLZ",
      "6VMI",
      "6XYW",
      "6ZM5",
      "6ZM6",
      "7A5F",
      "7A5G",
      "7A5H",
      "7A5I",
      "7A5J",
      "7A5K",
      "7L08",
      "7L20",
      "7O9K",
      "7O9M",
      "7ODR",
      "7ODS",
      "7ODT",
      "7OF0",
      "7OF2",
      "7OF3",
      "7OF4",
      "7OF5",
      "7OF6",
      "7OF7",
      "7OI6",
      "7OI7",
      "7OI8",
      "7OI9",
      "7OIA",
      "7OIB",
      "7OIC",
      "7OID",
      "7OIE",
      "7P2E",
      "7PD3",
      "7PKT",
      "7PNT",
      "7PNU",
      "7PNV",
      "7PNW",
      "7PNX",
      "7PNY",
      "7PNZ",
      "7PO0",
      "7PO1",
      "7PO2",
      "7PO3",
      "7PO4",
      "7QH6",
      "7QH7",
      "7QI4",
      "7QI5",
      "7QI6",
      "8A22",
      "8ANY",
      "8APN",
      "8APO",
      "8CSP",
      "8CSQ",
      "8CSR",
      "8CSS",
      "8CST",
      "8CSU",
      "8OIN",
      "8OIP",
      "8OIQ",
      "8OIR",
      "8OIS",
      "8OIT",
      "8PK0",
      "8QSJ"
    ]


# get trna cterm
def biopythin_chain_get_cterm_residues(bpchain):
    ...

# get mtrna PTC residues 
def res_nbhd(cterm_residues)->List[Residue]:
    ...

# get mtrna PTC residues, a good 10-15 of them -- whichver radius that works out to.
def res_nbhd(cterm_residues)->List[Residue]:
    ...

# get mtrna PTC residues, a good 10-15 of them -- whichver radius that works out to.
def project_residues(cterm_residues)->List[Residue]:
    ...


# T is a landmark with method project_into, project_from, data D and flag `present`
# basically assgin to every node of the taxonomy tree the the landmark with the data where there is one
# "project" from extant nodes to the rest preferring proximal nodes as sources
class GlobalTaxonomy[T]():
   
    ...


def get_rcsb_ids():
    return rcsb_ids