#!/bin/bash
# Re-runnable: set up ~/venvs/salterns_modelling on TF with the same ModelSEEDpy commit and the
# ModelSEEDDatabase files ModelSEEDpy's from_local() reads, both from git archives of the local checkouts
# cobra 0.29.1 pins pandas<3 but the local py_venv runs it with pandas 3.0.3, so cobra is installed --no-deps
# to match the local environment.
# (read-only on Andrew's repos: only `git archive` / `git rev-parse` are run).
set -euo pipefail
MSPY=$HOME/Documents/ModelSEEDpy
MSDB=$HOME/Documents/ModelSEEDDatabase
PY_COMMIT=$(git -C $MSPY rev-parse HEAD)          # sync-cshenry-main
DB_COMMIT=$(git -C $MSDB rev-parse HEAD)          # tautomer-adenine-n1_20261113
RUN=tf_runs/salterns_modelling
ssh TF "mkdir -p ~/$RUN/src/ModelSEEDpy ~/$RUN/src/ModelSEEDDatabase"
git -C $MSPY archive $PY_COMMIT | ssh TF "tar x -C ~/$RUN/src/ModelSEEDpy"
DBFILES=$(cd $MSDB && git ls-tree --name-only $DB_COMMIT Biochemistry/ | grep -E '^Biochemistry/(compound|reaction)_[0-9]+\.json$' | tr '\n' ' ')
git -C $MSDB archive $DB_COMMIT $DBFILES \
  Biochemistry/Aliases/Unique_ModelSEED_Compound_Aliases.txt Biochemistry/Aliases/Unique_ModelSEED_Compound_Names.txt \
  Biochemistry/Aliases/Unique_ModelSEED_Reaction_Aliases.txt Biochemistry/Aliases/Unique_ModelSEED_Reaction_Names.txt \
  Biochemistry/Aliases/Unique_ModelSEED_Reaction_ECs.txt Biochemistry/Structures/Unique_ModelSEED_Structures.txt \
  | ssh TF "tar x -C ~/$RUN/src/ModelSEEDDatabase"
printf "ModelSEEDpy\t%s\t%s\nModelSEEDDatabase\t%s\t%s\n" "$PY_COMMIT" "$(git -C $MSPY branch --show-current)" \
  "$DB_COMMIT" "$(git -C $MSDB branch --show-current)" | ssh TF "cat > ~/$RUN/src/COMMITS.tsv"
ssh TF "set -e; cd ~/$RUN; \
  [ -x ~/venvs/salterns_modelling/bin/python ] || ~/.local/bin/uv venv ~/venvs/salterns_modelling --python 3.12; \
  ~/.local/bin/uv pip install --python ~/venvs/salterns_modelling/bin/python -q \
    optlang==1.9.1 gurobipy==13.0.1 swiglpk==5.0.13 numpy==2.4.6 pandas==3.0.3 scipy==1.17.1 \
    sympy==1.14.0 networkx==3.6.1 chemicals==1.3.3 chemw==0.3.5 scikit-learn==1.7.0 matplotlib==3.10.3 \
    deepdiff==6.7.1 jinja2==3.1.6 requests==2.32.4 python-libsbml==5.20.5 depinfo==2.2.0 diskcache==5.6.3 \
    pydantic==2.11.9 rich==14.1.0 ruamel-yaml==0.18.15 appdirs==1.4.4 importlib-resources==6.5.2 httpx==0.28.1 \
    future==1.0.0 multiprocess==0.70.18 tqdm==4.67.1; \
  ~/.local/bin/uv pip install --python ~/venvs/salterns_modelling/bin/python -q --no-deps cobra==0.29.1 ./src/ModelSEEDpy; \
  ~/.local/bin/uv pip freeze --python ~/venvs/salterns_modelling/bin/python > src/pip_freeze.txt"
