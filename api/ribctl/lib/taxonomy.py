from ete3 import NCBITaxa

BACTERIA = 2
ARCHAEA  = 2157
EUKARYA  = 2759

ncbi = NCBITaxa()

def __to_names(d:dict):
    nu = {}
    for k in d.keys():
        if type(d[k]) == dict:
            nu[ [*ncbi.get_taxid_translator([int(k)]).values() ] [0]]=__to_names(d[k])
        else:
            nu[[* ncbi.get_taxid_translator([int(k)]).values() ][0]] = d[k]
    return nu





def __classify_profile(d:dict)->dict:
    s={}
    for _ in d['source_organisms']:
        if _ in s:
            s[_]+=1
        else:
            s[_]= 1

    h={}
    for _ in d['host_organisms']:
        if _ in h:
            h[_]+=1
        else:
            h[_]= 1

    top_org = [* sorted(s.items(), key=lambda x: x[1], reverse=True) ][0][0]

    return {
        "rcsb_id"       : d["rcsb_id"],
        "classified_as" : top_org,
        "src_orgs"      : s,
        "host_orgs"     : h,
    }