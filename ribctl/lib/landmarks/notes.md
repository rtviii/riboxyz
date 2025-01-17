## ptc_via_trna:

- for each ribosome type a structure is picked as a reference

- the trna residues at closest to ptc are pre-saved in `assets` 

- for a given target structure, we grab the landmarks on the references chains + _target_ rRNA 

- alignment by proxy of SequenceMappingContainer

- badabing badaboom, you got residues in target.