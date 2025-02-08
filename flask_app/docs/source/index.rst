.. Reptiles documentation master file, created by
   sphinx-quickstart on Tue Jun 11 13:53:12 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Reptiles's documentation!
====================================

Reptiles is a web app for automation of assembly design written in Python. It allows users to construct and design assemblies of 
multiple large fragments at once. The web app provides a friendly user interface and follows the Golden Braid rules for multipart 
iterative assemblies based on type 2s restriction enzymes and a standard grammar for easy part exchange and integration. The user 
uploads the sequences of the domesticated parts to be assembled in a genbank file format and then determines the size of the final 
construct in addition to the number of constructs to be designed. (In a case of multiple final constructs, the algorithm will design 
overlaps between them for further yeast assembly). The algorithm then creates an assembly plan for each construct and simulates the 
assemblies. As a result, the user receives assembly reports with graphs and additional information (warnings and summary reports). 
Moreover, a sequence genbank file for each construct and intermediate fragments is generated, which can be later visualized in 
Benchling or AddGene. In a few minutes, one can design and simulate several constructs of hundreds of kilobase pairs in size with 
extensive information on every intermediate step and sequence files. Reptiles aims to facilitate large-scale production of synthetic 
parts by reducing the time required for designing assembly plans and making it less error-prone.

.. toctree::
   usage
   :maxdepth: 2
   :caption: Contents:  
   