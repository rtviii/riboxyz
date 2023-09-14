
From today, longer and shorter term angles.

# "Expanding" the nomenclature

This nomenclature crusade should, ideally, driven by a concrete application, but it would be nice to at least have an account of useful/desirable categories and progress on each in one place. This exists in the form of a [ google sheet ](https://docs.google.com/spreadsheets/d/1mapshbn1ArofPN-Omu8GG5QdcwlJ0ym0BlN252kkUBU/edit?usp=sharing) at the moment. 

Obvious, but better stated explicitly: 

The ideal case for a subcomponent category is when it can be identified cheaply and reliably. The general utility of having these categories concretized and equipped with an unequivocal identification criterion is being able to parse them out from any structure with all the usual benefits of having a heterogeneity dataset like that. 


Some general issues with methodically parsing members of these categories are:

1. a lot of the criteria identifying them seem purely experimental and would be ad-hoc if translated directly into code.
2. the more granular categories just dissolve into bleeding edge of science and are still in flux, lack primary data etc.
 
With that said, there is definitely a subset of thes categories/types that has both scientific consensus and identifiable features. Basically, I see a balancing act between recruiting as many easily-identifiable types to the structural profile of the ribosome as opposed to gettig caught up in categorizing something ephemeral. The moment discrete categories start to fail and probably a few moments before that -- it's time to take look at shapes and density maps. I think this holds for both macro and micro features (ex. whole-struct assembly stages vs rProtein domain motions during elongation).


The categories so far (in no order of importance):

## Structural


- rRNAs
- rProteins
- tRNAs
- fProteins (Factors)
- Ligands

This things vary between species, of course. rProteins are covered by Ban and are basically "solved". rRNA likewise. Not clear how much classification tRNAs need, but a distinction should be made betwee fMET-type and AA-carrying trnas. Factors are a whole story onto themselves and are probably best addressed by lifecycle stage they belong to. Ligands are an "everything else" category, but medically-relevant targets can be chaffed out of it by comparing them against DrugBank (we have that).

2021 xu et al., elongation factors: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8807479/pdf/fmolb-08-816398.pdf

2018 frank and hashem, initiation factors: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6318078/ (edited) 

## Landmarks

- PTC
- Exit Tunnel Port
- SRL
- ASL
- *general rRNA Domains and Helicies*

Some of this stuff we get can probably get for free when banging out one of the "Structural" categories and others would be impossible to identify in the absence of some consensus sequence. Some of these landmarks will probably suffer from poorly resolved.


## Conformation/Lifecycle Stages

- Biogenesis
- Assembly
- Recycling

- Initiation
- Elongation
- Translocation
- Termination
- Release

Of these -- initiation, elongation, translocatio and termination are probably the best observed in the literature and can be somewhat lumped into one big "middle piece".
Assembly is generally complicated and is probably easier addressed through densitymap-based methods for now, but extensive reviews exist of its stages and the actors involved. Biogenesis is pretty much a nominal category for now, it seems: the structure is extremely nascent rRNA and we are treading on genetists' soil by this point; yet if at any point in the future we wanted to build bridges between structure and gene expression, this'll be the place to start.

I'll get a repo of going that keeps the literature observing each of the above.
https://github.com/rtviii/ribosome-resources

# General research

- Talk to Aryan re stages of elongation, work out the criteria for initiation, elongation, termination.

# Density Maps:

- Get all EMBD maps (start with E.coli) to start playing with them (rough classifcation/wasserstein distance). We can see about how all our tools can be brought to bear on assembly hence.

- See about plugging ChimeraX API into ribxz.

