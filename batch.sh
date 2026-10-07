#!/bin/bash
# Launch JupyterLab (dolfinx 0.7.2) with the repository mounted; open ch-limit-CG2.ipynb
# (figures of the paper, equal-order P2 elements) or ch-limit-RTDG.ipynb (Raviart-Thomas / DG pair).

docker run --rm --init -p 8888:8888 -v "$(pwd)":/root/shared dolfinx/lab:v0.7.2 ./shared/ch-limit-CG2.ipynb
