#!/bin/bash

python ./1.load.py

python ./4.af.py

python ./5_.generate_vep_variant_annots.py

python ./5_.add_db_vep_variant_annots.py

python ./7.nHetHom.py

python ./8.drop_samples.py

python ./9.export.ssvs.py

python ./9.export.mapd.py
