#!/bin/bash -e

awk '$5 >= -0.1 && $5 <= 0.1'  controls_calc_z.tsv > controls_z_near_zero.tsv
