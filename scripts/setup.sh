#!/usr/bin/env bash

python3 -m pybel_tools manage role add "admin"
python3 -m pybel_tools manage role add "scai"

python3 -m pybel_tools manage user add --admin "cthoyt@gmail.com" "pybeladmin"
python3 -m pybel_tools manage user add --scai "scai@example.com" "pybeltest"
python3 -m pybel_tools manage user add "test@example.com" "pybeltest"
