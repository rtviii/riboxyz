import argparse
import prody as prd
from prody import calcShannonEntropy

# parser = argparse.ArgumentParser (description='Ribosome Cli')

# #TODO:
# # - SEPARATE LIB CODE from CLI CODE
# # - retype the typescript profile code in python

# parser.add_argument ("-it"       , "--itern"               , type= int      ,                 help = "The number of iterations"                                                                                            )
# # parser.add_argument ("-itstart"  , "--iter_start"              , type   = int          ,required =True,                 help = "The number of iterations"                                                                                                                                )
# # parser.add_argument ("-itend"    , "--iter_end"                , type   = int          ,required =True,                 help = "The number of iterations"                                                                                                                                )
# args                         = parser .parse_args()

# GENERATION 					 = 1000
# ITSTART                      = int(args.iter_start)
# ITEND                        = int(args.iter_end)
# REPLICATE_N                  = int (args .siminst if args.siminst is not None else 0)
# OUTDIR                       = args.outdir if args.outdir is not None else 0
# RESSURECT_PATH               = args.resurrect if args.resurrect is not None else 0



def barr2str (bArr):

    return ''.join([ x.decode("utf-8") for x in bArr])
    

# open msa file with prody
msa_file = prd.parseMSA('PV_uL2.fasta')
s        = calcShannonEntropy(msa_file, omitgaps=True)
print(sum(s))





