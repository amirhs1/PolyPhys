#!/bin/bash

lmp_serial -log minimize_dna.log < minimize_dna.lmp
lmp_serial < input.lmp