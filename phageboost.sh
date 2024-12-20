#!/bin/bash

base=$(basename "$0")
dir=$(dirname "$0")
path="$dir/PhageBoost"
export PYTHONPATH="$path:$PYTHONPATH"

