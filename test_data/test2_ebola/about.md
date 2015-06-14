#Ebola Search

This test is from a BLAST search I performed on Ebola data before and after/during the west african epidemic 2014. The query data is from viprbrc.org (which gets the data from NCBI) of an 1976 ebola case in Yambuku, DRK.

First search was performed with MegaBLAST on the dutch [Andromeda Server](http://galaxy.nbic.nl/) with an NCBI nuleotide database from 21.05.2013, the second on the [Galaxy Server of the Georgia Aquarium Whale shark project](http://whaleshark.georgiaaquarium.org/) with a database from 20.10.2014.

The results were obtained from running the main.py script via the command line:

´´´bash
python main.py -o ebola_old.tabular -n ebola_new.tabular -Al -s result.txt
´´´