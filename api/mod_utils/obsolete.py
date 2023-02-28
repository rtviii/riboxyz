from operator import sub
from numpy import log
from rest_framework.decorators import api_view
from rest_framework.response import Response
from neo4j import  Result, GraphDatabase

from rest_framework.serializers import Serializer
import time
import os
#-⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅,⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯
uri        =  os.getenv( 'NEO4J_URI'                                      )
authglobal = (os.getenv( 'NEO4J_USER'      ),os.getenv( 'NEO4J_PASSWORD' ))
current_db =  os.getenv( 'NEO4J_CURRENTDB'                                )



@api_view(['GET'])
def get_struct(request):
    params = dict(request.GET)
    pdbid  = str.upper(params['pdbid'][0])
    cypher = """
    match (n:RibosomeStructure{{rcsb_id:"{pdbid}"}})
    optional match (rr:RNA)-[]-(n)
    with n, collect(rr) as rrna
    optional match (rp:Protein)-[]-(n)
    with n, rrna,  collect(rp) as rps
    optional match (l:Ligand)-[]-(n)
    with n, rrna, rps, collect(l) as ligs
    return {{structure: n, ligands: ligs,rnas: rrna, proteins: rps}}
    """.format_map({"pdbid":pdbid})
    return Response(_neoget(cypher))


@api_view(['GET' ])
def get_all_ligands(request): 
    CYPHER_STRING="""
        match (l:Ligand)-[]-(r:RibosomeStructure)  where 
        not l.chemicalName  contains "ION" 
        and not l.chemicalName contains "CLUSTER"
        and not l.chemicalName contains "["
        and r.expMethod <> "X-RAY DIFFRACTION"
        return {{  
        polymer    : false,
        description: l.chemicalName,
        chemicalId : l.chemicalId,
        presentIn  : {{
                src_organism_ids: r.src_organism_ids,
                description     : l.chemicalName,
                citation_title  : r.citation_title,
                expMethod       : r.expMethod,
                rcsb_id         : r.rcsb_id,
                resolution      : r.resolution
            }}
        }}
    """.format()
    return Response(_neoget(CYPHER_STRING))

@api_view(['GET' ])
def get_all_ligandlike(request):
    CYPHER_STRING = """
        match (l {{ligand_like:true}})-[]-(r:RibosomeStructure) 
        where r.expMethod <> "X-RAY DIFFRACTION"
        return {{
            polymer     : true,
            description : l.rcsb_pdbx_description,
            presentIn  : {{
                auth_asym_id    : l.auth_asym_id,
                src_organism_ids: r.src_organism_ids,
                description     : l.rcsb_pdbx_description,
                citation_title  : r.citation_title,
                expMethod       : r.expMethod,
                rcsb_id         : r.rcsb_id,
                resolution      : r.resolution
            }}
        }}""".format()
    return Response(_neoget(CYPHER_STRING))

#? ---------------------------STRUCTS

@api_view(['GET']) 
def get_ligands_by_struct(request):
    CYPHER_STRING="""match (n:RibosomeStructure)-[]-(l:Ligand)
           return {{ title: n.citation_title, struct:n.rcsb_id, organism:n.src_organism_names, taxid:n.src_organism_ids, 
           ligands:collect({{ chemid: l.chemicalId, name:l.chemicalName, number:number_of_instances }})}}""".format_map({})
    return Response(_neoget(CYPHER_STRING))


    
    
    

@api_view(['GET'])
def get_all_structs(request):
    start = time.time()
    print("Begun fetching all structures.")
    CYPHER_STRING="""match (ribs:RibosomeStructure) 
        unwind ribs as rb

        optional match (l:Ligand)-[]-(rb)
        with collect(l.chemicalId) as ligs, rb

        optional match (rps:Protein)-[]-(rb)
        with ligs, rb, collect({{
            auth_asym_id                   : rps.auth_asym_id,
            nomenclature                   : rps.nomenclature,
            entity_poly_seq_one_letter_code: rps.entity_poly_seq_one_letter_code
            }}) as rps

        optional match (rnas:RNA)-[]-(rb)
        with ligs, rb, rps, collect({{
            auth_asym_id                   : rnas.auth_asym_id,
            nomenclature                   : rnas.nomenclature,
            entity_poly_seq_one_letter_code: rnas.entity_poly_seq_one_letter_code
            }}) as struct_rnas

        return {{
            struct : rb         ,
            ligands: ligs       ,
            rps    : rps        ,
            rnas   : struct_rnas
            }} 
        """.format()
    qres = _neoget(CYPHER_STRING)
    print("Returning all structures. Query took {} s".format(time.time()-start) )
    return Response(qres)
#? ---------------------------PROTEINS

@api_view(['GET'])
def get_banclass_for_chain(request):
    params       = dict(request.GET)
    pdbid        = str.upper(params['pdbid'][0])
    auth_asym_id = str(params['auth_asym_id'][0])
    cypher       = """match (n:RibosomeStructure {{rcsb_id:"{pdbid}"}})-[]-(c:Protein{{auth_asym_id:"{auth_asym_id}"}})-[]-(pc:ProteinClass) return pc.class_id
    """.format_map({"pdbid":pdbid,"auth_asym_id":auth_asym_id})
    return Response(_neoget(cypher))

@api_view(['GET'])
def get_banclasses_metadata(request):
    params  = dict(request.GET)
    family  = str(params['family'][0]).lower() # u e b
    subunit = str(params['subunit'][0]).lower()

    if subunit == "ssu":
        fstring = 'toLower(n.class_id) contains "s" or toLower(n.class_id) contains "bthx" or toLower(n.class_id) contains "rack"' 
    elif subunit == "lsu": 
        fstring = 'toLower(n.class_id) contains "l"' 

    CYPHER_STRING="""
    match (n:ProteinClass)-[]-(rp:Protein)-[]-(s:RibosomeStructure) where  toLower(n.class_id) contains "{}"  and {} 
    unwind s.`src_organism_ids` as orgid
    with collect(distinct orgid) as allorgs, n as n, collect(s.rcsb_id) as structures, collect(distinct rp.pfam_comments) as comments
    return {{banClass: n.class_id, organisms:  allorgs, comments:comments, structs: structures }}""".format(family, fstring)

    return Response(_neoget(CYPHER_STRING))

@api_view(['GET'])
def list_nom_classes(request):
    CYPHER_STRING="""
    match (b:ProteinClass)-[]-(rp)-[]-(str:RibosomeStructure)
    with str, b, rp
    return {{
    nom_class: b.class_id,
    rps:collect({{
        organism_desc: rp.src_organism_names,
        organism_id  : rp.src_organism_ids,
        uniprot      : rp.uniprot_accession,
        parent       : str.rcsb_id,
        parent_reso  : str.resolution,
        strand_id    : rp.entity_poly_strand_id
        }}),
    presentIn:collect(str.rcsb_id)}}""".format_map({})
    return Response(_neoget(CYPHER_STRING))

@api_view(['GET'])
def gmo_nom_class(request):
    params = dict(request.GET)
    ban    = str(params['banName'][0])
    CYPHER_STRING="""
    match (rib:RibosomeStructure)-[]-(n:Protein)-[]-(nc:ProteinClass{{class_id:"{ban}"}})

    return {{  
    parent_resolution                  : rib.resolution,
    parent_year                        : rib.citation_year,
    parent_method                      : rib.expMethod,
    parent_rcsb_id                     : n.parent_rcsb_id,
    pfam_accessions                    : n.pfam_accessions,
    pfam_comments                      : n.pfam_comments,
    pfam_descriptions                  : n.pfam_descriptions,
    src_organism_ids                   : n.src_organism_ids,
    src_organism_names                 : n.src_organism_names,
    uniprot_accession                  : n.uniprot_accession,
    rcsb_pdbx_description              : n.rcsb_pdbx_description,
    entity_poly_strand_id              : n.entity_poly_strand_id,
    entity_poly_seq_one_letter_code    : n.entity_poly_seq_one_letter_code,
    entity_poly_seq_one_letter_code_can: n.entity_poly_seq_one_letter_code_can,
    entity_poly_seq_length             : n.entity_poly_seq_length,
    entity_poly_polymer_type           : n.entity_poly_polymer_type,
    entity_poly_entity_type            : n.entity_poly_entity_type,
    surface_ratio                      : n.surface_ratio,
    nomenclature                       : n.nomenclature
        }}
    """.format_map({
        "ban":ban
    })

    return Response(_neoget(CYPHER_STRING))

@api_view(['GET']) 
def banclass_annotation(request):
    params   = dict(request.GET)
    banclass = str(params['banclass'][0])

    CYPHER_STRING=f"""
            match (n:ProteinClass{{class_id:"{banclass}"}})-[]-(rp:Protein) 
            with rp.rcsb_pdbx_description as dd 
            return  dd limit 6;
           
           """
    return Response(_neoget(CYPHER_STRING))

@api_view(['GET'])
def nomclass_visualize(request):

    params = dict(request.GET)
    ban    = str(params['ban'][0])

    CYPHER_STRING="""
    match (n:ProteinClass)-[]-(rp:Protein) where n.class_id ="{}" return  {{
    class: n.class_id,
    members: collect({{parent: rp.parent_rcsb_id, chain:rp.entity_poly_strand_id}}),
    comments:collect(distinct rp.pfam_comments)}} """.format(ban)
    
    return Response(_neoget(CYPHER_STRING))

@api_view(['GET'])
def proteins_number(request):
    CYPHER_STRING="""match (n:Protein) return count(n);"""
    return Response(_neoget(CYPHER_STRING))


#? ---------------------------RNA
                     
@api_view(['GET']) 
def get_rnas_by_struct(request):
    CYPHER_STRING="""match (n:RibosomeStructure)-[]-(r:RNA) where toLower(r.rcsb_pdbx_description)  contains "mrna" or toLower(r.rcsb_pdbx_description) contains "trna"
        or toLower(r.rcsb_pdbx_description)  contains "m-rna" or toLower(r.rcsb_pdbx_description)  contains "t-rna"  or toLower(r.rcsb_pdbx_description)  contains "messenger" or toLower(r.rcsb_pdbx_description)  contains "transfer"
        return {{struct:n.rcsb_id, rnas:collect(r.rcsb_pdbx_description)}};""".format_map({})
    return Response(_neoget(CYPHER_STRING))

@api_view(['GET']) 
def get_rna_class(request):
    params        = dict(request.GET)
    rna_class     = str(params['rna_class'][0])
    CYPHER_STRING  = """
    match (c:RNAClass {{ class_id:"{}" }})-[]-(n)-[]-(rib:RibosomeStructure)
    return {{
        parent_year                         : rib.citation_year                      ,
        parent_resolution                   : rib.resolution                         ,
        parent_citation                     : rib.citation_title                     ,
        parent_method                       : rib.expMethod                          ,
        asym_ids                            : n.asym_ids                           ,
        auth_asym_id                        : n.auth_asym_id                       ,
        nomenclature                        : c.class_id                           ,
        parent_rcsb_id                      : n.parent_rcsb_id                     ,
        src_organism_names                  : n.src_organism_names                 ,
        host_organism_names                 : n.host_organism_names                ,
        src_organism_ids                    : n.src_organism_ids                   ,
        host_organism_ids                   : n.host_organism_ids                  ,
        rcsb_pdbx_description               : n.rcsb_pdbx_description              ,
        entity_poly_strand_id               : n.entity_poly_strand_id              ,
        entity_poly_seq_one_letter_code     : n.entity_poly_seq_one_letter_code    ,
        entity_poly_seq_one_letter_code_can : n.entity_poly_seq_one_letter_code_can,
        entity_poly_seq_length              : n.entity_poly_seq_length             ,
        entity_poly_polymer_type            : n.entity_poly_polymer_type           ,
        entity_poly_entity_type             : n.entity_poly_entity_type            ,
        ligand_like                         : n.ligand_like                        
    }}
    """.format(rna_class)

    return Response(_neoget(CYPHER_STRING))

@api_view(['GET']) 
def ranged_align(request):
    params = dict(request.GET)
    print("-------------------+------------------")
    print("GOT PARAMS", params)
    print("-------------------+------------------")

    r1start = int(params['r1start'][0])
    r1end   = int(params['r1end'][0])

    r2start = int(params['r2start'][0])
    r2end   = int(params['r2end'][0])

    struct1       = params['struct1'][0].upper()
    struct2       = params['struct2'][0].upper()
    auth_asym_id1 = params['auth_asym_id1'][0]
    auth_asym_id2 = params['auth_asym_id2'][0]

    RANGED_ALIGNMENT_SCRIPT= os.path.join(str( PROJECT_PATH ), 'static_files','ranged_align.py')
    os.system("python3 {} {} {} {} {} {}-{} {}-{}".format(RANGED_ALIGNMENT_SCRIPT,struct1,struct2, auth_asym_id1, auth_asym_id2, r1start,r1end, r2start,r2end))

    alignedfile = os.environ["TEMP_CHAIN"]

    try:
        doc = open(alignedfile, 'rb')
    except: 
        print(f"Could not find {alignedfile}. Exited")
        return Response(-1)

    response = HttpResponse(FileWrapper(doc), content_type='chemical/x-mmcif')
    response['Content-Disposition'] = 'attachment; filename="{}-{}_{}-{}.cif"'.format(struct1,auth_asym_id1,struct2,auth_asym_id2)
    return response

@api_view(['GET',])
def get_chain(request):
    params   = dict(request.GET)
    chainid  = params['chainid'][0]
    structid = str.upper(params['structid'][0])
    filename = "{}_subchain_{}.pdb".format(structid, chainid)

    file_handle = os.path.join(STATIC_ROOT,structid, filename)

    document = open(file_handle, 'rb')
    response = HttpResponse(FileWrapper(document), content_type='chemical/x-pdb')
    response['Content-Disposition'] = 'attachment; filename="{}_subchain_{}.pdb"'.format(structid, chainid)
    return response

@api_view(['GET',])
def download_ligand_nbhd(request):
    params   = dict(request.GET)
    structid = params['structid'][0].upper()
    chemid   = params['chemid'][0].upper()

    filename   = "LIGAND_{}.json".format(chemid)
    filehandle = os.path.join(STATIC_ROOT, structid, filename)

    try:
        doc = open(filehandle, 'rb')
    except: 
        print(f"Could find {filehandle}.Exited")
        return Response(-1)

    response = HttpResponse(FileWrapper(doc), content_type='application/json')
    response['Content-Disposition'] = 'attachment; filename="{}_LIGAND_{}.json"'.format(structid, chemid)
    return response

@api_view(['GET',])
def ligand_prediction(request):

    params     = dict(request.GET)
    # it is either a chemical id (if is_polymer == False) 
    # or a entity_poly_strand_id in the case of a ligand-like polymer (is_polymer  ==True)
    ligandlike_id = params['ligandlike_id'][0]
    src_struct    = params['src_struct' ][0].upper()
    tgt_struct    = params['tgt_struct' ][0].upper()
    is_polymer    = str( params['is_polymer' ][0] )
    print("Attempting to render  {} from {}(orig) in {}.".format(ligandlike_id,src_struct,tgt_struct))
    prediction_filename = "PREDICTION_{}_{}_{}.json".format     (ligandlike_id,src_struct ,tgt_struct          )
    filehandle          = os.path  .join(STATIC_ROOT, tgt_struct, prediction_filename)

    orig_binding_site_handle = os.path.join(STATIC_ROOT, src_struct, "{}_{}.json".format("POLYMER" if is_polymer.lower() == "true" else "LIGAND", ligandlike_id))
    target_json_handle       = os.path.join(STATIC_ROOT, tgt_struct, "{}.json".format(tgt_struct))
    with open(orig_binding_site_handle, 'rb') as infile:
        bsite = BindingSite(json.load(infile))

    with open(target_json_handle, 'rb') as target_handle:
        target_handle   = json.load(target_handle)

    #* Transpose Ligand Script
    prediction = transpose_ligand.init_transpose_ligand(src_struct,tgt_struct, target_handle, bsite)
    return Response(prediction)

@api_view(['GET',])
def get_ligand_nbhd(request):
    params        = dict(request.GET)
    print("----------------PARAMS ")
    print(params)
    src_struct    = params['src_struct'][0].upper()
    ligandlike_id = params['ligandlike_id'][0]
    is_polymer    = str( params['is_polymer'][0] )
    filehandle    = os.path .join(STATIC_ROOT, src_struct, "{}_{}.json".format("POLYMER" if is_polymer.lower() == 'true' else "LIGAND", ligandlike_id))

    print(f"Returning data from file {filehandle}")
    try:
        with open(filehandle, 'rb') as infile:
            data = json.load(infile)
            return Response(data)
    except error: 
        print("errored out", error)
        
        return Response(-1)

@api_view(['GET', ])
def cif_chain(request):
    params     = dict(request.GET)
    chainid    = params['chainid'][0].upper()
    struct     = params['structid'][0].upper()

    filename   = "{}_STRAND_{}.cif".format(struct,chainid)
    filehandle = os.path.join(STATIC_ROOT, struct,'CHAINS', filename)

    try:
        doc = open(filehandle, 'rb')
    except: 
        return Response("File not found")

    response = HttpResponse(FileWrapper(doc), content_type='chemical/x-mmcif')
    response['Content-Disposition'] = 'attachment; filename="{}"'.format(filename)

    return response


@api_view(['GET'])
def download_structure(request):
    params     = dict(request.GET)
    struct_id    = params['struct_id'][0].upper()

    filename   = "{}_modified.cif".format(struct_id)
    filehandle = os.path.join(STATIC_ROOT, struct_id, filename)

    try:
        doc = open(filehandle, 'rb')
    except: 
        return Response("File not found")
    response = HttpResponse(FileWrapper(doc), content_type='chemical/x-mmcif')
    response['Content-Disposition'] = 'attachment; filename="{}"'.format(filename)
    return response

@api_view(['GET', ])
def cif_chain_by_class(request):
    params     = dict(request.GET)
    classid    = params['classid'][0]
    struct     = params['struct'][0].upper()

    CYPHER = """match (n:RibosomeStructure)-[]-(r)-[]-(b) where n.rcsb_id ="{}" and b.class_id = "{}"
    return {{ struct: n.rcsb_id, auth_asym_id: r.auth_asym_id }}""".format(struct,classid)

    chains = _neoget(CYPHER)
    if len( chains ) < 1 :
        return Response("Not found")
    auth_asym_id     = chains[0]['auth_asym_id']
    filename   = "{}_STRAND_{}.cif".format(struct,auth_asym_id)
    filehandle = os.path.join(STATIC_ROOT, struct.upper(), "CHAINS", filename)

    try:
        doc = open(filehandle, 'rb')
    except: 
        return Response("File not found")

    response = HttpResponse(FileWrapper(doc), content_type='chemical/x-mmcif')
    response['Content-Disposition'] = 'attachment; filename="{}"'.format(filename)
    return response


def tunnel(request):
    params     = dict(request.GET)
    struct     = params['struct'][0].upper()
    filetype   = params['filetype'][0]

    cetrline_filehandle  =  os.path.join(STATIC_ROOT, struct,'TUNNEL', 'csv', 'centerline.csv')
    report_filehandle    =  os.path.join(STATIC_ROOT, struct,'TUNNEL', f"{struct}_TUNNEL_REPORT.json")


    if filetype== 'report':
        print("GOT REQUEST FOR REPORT with params", params)
        try:
            doc = open(report_filehandle, 'rb')
        except: 
            return Response("File not found")

        response = HttpResponse(FileWrapper(doc), content_type='application/json')
        response['Content-Disposition'] = 'attachment; filename="{}"'.format(f"{struct}_tunnel_report.json")

        return response
    elif filetype == 'centerline':
        try:
            doc = open(cetrline_filehandle, 'rb')
        except: 
            return Response("File not found")

        response = HttpResponse(FileWrapper(doc), content_type='text/csv')
        response['Content-Disposition'] = 'attachment; filename="{}"'.format(f"{struct}_tunnel_centerline.csv")

        return response

@api_view(['GET', ])
def downloadArchive(request):

    params  = dict(request.GET)
    print("received dict " , params)

    if 'rna' in params:
        rnas    = params['rna']
        rnas  =[*map(lambda x: x.split('.'),rnas)]
        rnas  =[*map(lambda x:  os.path.join(STATIC_ROOT,x[0].upper(),'CHAINS',"{}_STRAND_{}.cif".format(x[0].upper(),x[1])),rnas)]
        print("dict:", rnas)
    else:
        rnas = []

    if 'structs' in params:
        structs = params['structs']
        structs = [*map(lambda x:  os.path.join(STATIC_ROOT,x.upper(),"{}.cif".format(x.upper())),structs)]
        print("structs", structs)
    else:
        structs = []

    file_names = [*structs,*rnas]
    zip_subdir = 'ribosome_xyz_archive'
    zf         = zipfile.ZipFile('temp_zip.zip', "w")

    for fpath in file_names:
        fdir, fname = os.path.split(fpath)
        zip_path = os.path.join(zip_subdir, fname)
        try:
            zf.write(fpath, zip_path)
            print(f"Failed to find file { fpath }")
        except:
            continue

    zf.close()
    r = open('temp_zip.zip','rb')
    os.remove('temp_zip.zip')
    return FileResponse(r)



# # TODO: nomenclature 
# @api_view(['GET'])
# def nomenclature(request):
#     params        = dict(request.GET)
#     if 'rcsb_id' in params and len( params['rcsb_id'] ) > 0:
#         cypher_single ="""match (n:RibosomeStructure{{rcsb_id:"{}"}})-[]-(c) where c:Protein or c:RNA 
#         return {{struct:n.rcsb_id, strand:c.auth_asym_id,nomenclature:c.nomenclature}}""".format(params['rcsb_id'][0].upper())
#         maps=  _neoget(cypher_single)
#         d={}
#         for _ in maps:
#             d.update({ _['auth_asym_id']:_['nomenclature'] })
#         return Response(d)
#     else:
#         cypher_all ="""match (n:RibosomeStructure)-[]-(c) where c:Protein or c:RNA 
#         return {{struct:n.rcsb_id, auth_asym_id:c.auth_asym_id,nomenclature:c.nomenclature}}""".format()
#         all_maps = {

#         }

#         resp = _neoget(cypher_all)
#         for _ in resp:
#             if _['struct'] in all_maps:
#                 all_maps[_['struct']].update({_['auth_asym_id']:_['nomenclature']})
#             else:
#                 all_maps.update(
#                     {
#                     _['struct']:{_['auth_asym_id']:_['nomenclature']}
#                     }
#                 )

#         return Response(all_maps)
