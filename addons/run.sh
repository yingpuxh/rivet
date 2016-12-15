#!bin/bash

mkfifo fifo.hepmc
less example.hepmc.gz  >> fifo.hepmc & 
rivet -a FatHiggsTagging --process higgs -o example.yoda fifo.hepmc


