#!/usr/bin/env python3
"""Thin shim – delegates to the ``useq2tracks`` package.

Nextflow automatically adds its workflow ``bin/`` dir to ``$PATH``,
so this script is resolved by name inside NF processes.  It requires
``useq2tracks`` to be installed in the execution environment
(``pip install useq2tracks``).
"""
from useq2tracks.check_samplesheet import main

if __name__ == "__main__":
    main()
