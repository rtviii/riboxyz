import typing
import requests
from  ribctl.lib.types.types_polymer import  list_LSU_Proteins,list_SSU_Proteins
import sys


B            = "2"
E            = "2759"
A            = "2157"
PROTEOVISION = lambda SU,_class, taxids: "https://ribovision3.chemistry.gatech.edu/showAlignment/{}/{}/{}".format(SU,_class,taxids)

def save_XSU_alns(SU:typing.Literal["SSU","LSU"]):
    protnames = list_SSU_Proteins if SU == "SSU" else list_LSU_Proteins 
    for protein in protnames:
        spec_letters, class_digit = protein.split("S" if SU == "SSU" else "L")
        if int(class_digit) < 10:
            class_digit = "0{}".format(class_digit)

        def tax_string(spec_letters):
            tax_string = ""
            if 'e' in spec_letters:
                return E
            if 'b' in spec_letters:
                return B
            if 'u' in spec_letters:
                return ",".join([B,E,A])
            return tax_string[:-1]

        URI         = PROTEOVISION(SU,spec_letters + ( "S" if SU == "SSU" else "L" ) + class_digit, tax_string(spec_letters))
        resp        = requests.get(URI)
        start,end   = resp.text.find("<pre>"), resp.text.find("</pre>")
        fasta_lines = resp.text[start+len("<pre>"):end].replace("&gt;",">")
        filename    = "{}_ribovision.fasta".format(protein)

        print("{}\t| {} \t| URI: {}".format(protein,resp.status_code, URI) )

        if resp.status_code == 200:
            with open(filename,"w") as f:
                f.write(fasta_lines)
                print("Wrote {} to {}".format(protein,filename))
         
save_XSU_alns("SSU")