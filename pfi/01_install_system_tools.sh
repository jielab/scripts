#!/usr/bin/env bash


# 🚩 01_install_system_tools
set -euo pipefail

sudo apt-get update
sudo apt-get install -y wget curl aria2 unzip gzip tar libxml2-utils build-essential r-base r-base-dev
