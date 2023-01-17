var fs = require('fs');

// I'm opening the file on the filesystem and parsing it as JSON
var source = JSON.parse(fs.readFileSync('./4v6u_67_B1.out.json', 'utf8'));

var RCSB_ID  = '4v6u';  // rcsb_id
var model    = '1';     // i'm assuming this is the model number 
var chain_id = "B1";    // chain id


interface rna_nucleotide {
    resnum: number,
    path  : number[]
}

const arr2str = (arr: number[]) => "M" + arr.join(',');
const svg_paths = source[RCSB_ID][model][chain_id]['rna_nucleotides']
    .reduce((acc: string[][], next: rna_nucleotide) => {
        return [...acc, arr2str(next.path)]
    }, [])