#!/bin/bash
set -euo pipefail

python ./5_.generate_vep_variant_annots.py

python ./5_.add_db_vep_variant_annots.py
