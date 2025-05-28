#!/bin/bash

source .venv/bin/activate
uv pip uninstall $PWD
git pull
uv pip install -e .
exit 0
