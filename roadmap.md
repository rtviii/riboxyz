


<!-- TODO: -->
- split ligand classes into instances and classes


# Features to stability:

- api
- url fixes

# DevOps:

- one-click docker deployment 
- how can we benefit from github actions?

# Later:

- "last updated" section
- automatic update per struct -- cron job every 24 weeks

- [x] rewrite struct init scripts in python
- [x] rewrtie cli in python
- [ ] graph.ribosome.xyz
- monitoring

# LIGANDS

- ligands/ligandlike should have a common category in the db. filter the ions too. some ligands are not redndered at all (ex. 5AFI.FME)

- gene ontology + bird + CCD

## Branching out to other resources
- interpro

## Static Files Server
- https://docs.djangoproject.com/en/4.1/howto/static-files/deployment/


## Integrate the PTC/conserved sites scripts scripts into the processing pipeline



# -------------- nlp stuff


Rasa or Dialogflow

 te a BERT or SciBERT model on a set of ribosome-related texts, and then use the model to extract information such as the names of subcomponents of the ribosome, the sequences of related genes, or the structures of ligands or antibiotics that interact with the ribosome. You could then use this information to generate code snippets that perform tasks like querying a database of ribosome-related information or visualizing the structure of the ribosom

 paCy for NER, NetworkX for graph analysis, and scikit-learn for clustering.



 Sure! To create a graph representation of the different components of the ribosome, you could start by using NLP techniques to extract relevant information from the text. Here are some steps you could follow:

    Preprocessing: Clean and preprocess the text data by removing stop words, punctuation, and other non-relevant information. You can use Python libraries like NLTK or SpaCy for this.

    Named Entity Recognition: Use a named entity recognition (NER) model to extract mentions of different ribosome components from the text. You can train your own NER model using libraries like spaCy or use pre-trained models like the BioBERT model, which has been specifically trained on biomedical texts.

    Dependency Parsing: Use a dependency parser to extract the relationships between the different ribosome components. The parser will identify the grammatical relationships between the different words in a sentence, which can help you identify which components are related to each other. You can use libraries like Stanford CoreNLP or spaCy for this.

    Graph Creation: Once you have extracted the ribosome components and their relationships from the text, you can create a graph representation of the components and their relationships. You can use a graph library like NetworkX or Neo4j for this.

Here's an example of how you might use these techniques to create a graph representation of the ribosome components:

    Preprocess the text data by removing stop words and other non-relevant information.

    Use a pre-trained NER model like BioBERT to extract mentions of different ribosome components from the text. For example, you might extract mentions of "small ribosomal subunit" or "rRNA".

    Use a dependency parser to extract the relationships between the different ribosome components. For example, you might identify that the "small ribosomal subunit" is a component of the "ribosome" and that it interacts with other components like "mRNA" and "tRNA".

    Create a graph representation of the ribosome components and their relationships using a library like NetworkX or Neo4j. Each component can be represented as a node in the graph, and the relationships between the components can be represented as edges.

Once you have created the graph representation of the ribosome components, you can use graph analysis techniques like clustering and community detection to group related components together and identify key substructures within the graph.

I hope this helps! Let me know if you have any other questions.