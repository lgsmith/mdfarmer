from __future__ import annotations

import json
import re
import subprocess as sp
import time
from os import environ
from pathlib import Path

import openmm as mm
import openmm.app as app
import openmm.unit as u
import parmed
