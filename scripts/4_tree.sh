#!/bin/bash

mkdir Tree


grapetree -p reads/SRS/Amend/*phylo_w12.plainIDs.fasta > Tree/SRS_tree.nwk

grapetree -p reads/LRS/Amend/*phylo_w12.plainIDs.fasta > Tree/LRS_tree.nwk

grapetree -p reads/Hybrid/Amend/*phylo_w12.plainIDs.fasta > Tree/Hybrid_tree.nwk


