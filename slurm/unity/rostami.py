#!/usr/bin/env python3
"""Add Rostami campaign provenance requirements to the immutable stage controller."""
import sys
import stabilization

stabilization.REQUIRED_LAUNCH_FILES |= {
    'scripts/make_rostami.jl', 'slurm/unity/rostami.sh', 'slurm/unity/rostami.py',
}

if __name__ == '__main__':
    try:
        stabilization.main()
    except (ValueError, OSError, KeyError) as error:
        print('saturn Rostami: ' + str(error), file=sys.stderr)
        sys.exit(1)
