import argparse
import os
import subprocess
import sys
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


def muscle_combine_profiles(msa_path1: str, msa_path2: str, out_filepath: str):
    """Combine two MSA-profiles into a single one. Used here to "append" a target sequence two the ribovision alignment. """
    cmd = ['/home/rxz/dev/docker_ribxz/cli/scripts/muscle3.8', '-profile','-in1', msa_path1, '-in2', msa_path2, '-out', out_filepath]
    subprocess.Popen(cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE, env=os.environ.copy()).wait()
    sys.stdout.flush()

def barr2str (bArr):
    return ''.join([ x.decode("utf-8") for x in bArr])
    

# open msa file with prody
msa_file = prd.parseMSA('PV_uL2.fasta')
s        = calcShannonEntropy(msa_file, omitgaps=False)
print("Unpetrurbed entropy: {}".format(sum(s)))
suspects = [ 
'bL34',
'bS2',
'uL2',
'uS3',
'uS5'
 ]

paths = [*map(lambda _: "5afi_unclassified_{}.fasta".format(_),suspects)]
for control in suspects:
    control_path = (lambda _: "5afi_unclassified_{}.fasta".format(_))(control)
    muscle_combine_profiles('PV_uL2.fasta', control_path, 'PV_uL2_with_{}.fasta'.format(control))


def compare_entropies(msa_path1: str, msa_path2: str):
    msa_file = prd.parseMSA('PV_uL2.fasta')
    msa_file_ctl = prd.parseMSA('PV_uL2.fasta')