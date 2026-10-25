#!/usr/bin/env python3
"""Check public static assets without modifying files, Git, or remote services.

This catches accidental model dumps, links, oversized assets, and common secret
formats. It is not a replacement for reviewing everything intended for release.
"""

import argparse
from pathlib import Path
import re
import stat
import sys


FILE_LIMIT = 90_000_000
TOTAL_LIMIT = 1_000_000_000
EXTENSIONS = {'.html', '.css', '.js', '.json', '.svg', '.png', '.jpg', '.jpeg',
              '.gif', '.webp', '.mp4', '.webm', '.csv', '.txt', '.md', '.ico', '.toml'}
COMPLETION_MARKER = b'Static gallery build completed.\n'
SECRET_NAMES = {'credentials', 'secrets', 'id_rsa', 'id_ed25519', 'id_ecdsa'}
SECRET_PATTERN = re.compile(
    rb'-----BEGIN (?:[A-Z0-9 ]+ )?PRIVATE KEY(?: BLOCK)?-----'
    rb'|(?:gh[pousr]_[A-Za-z0-9]{30}|github_pat_[A-Za-z0-9_]{30})'
    rb'|AKIA[0-9A-Z]{16}'
)


def check_gallery(root):
    root = Path(root)
    if root.is_symlink() or not root.is_dir():
        raise ValueError('gallery must be an ordinary directory')
    if not (root / 'index.html').is_file():
        raise ValueError('gallery/index.html is missing')
    count, total = 0, 0
    for path in sorted(root.rglob('*')):
        relative = path.relative_to(root)
        info = path.lstat()
        if stat.S_ISLNK(info.st_mode):
            raise ValueError('symbolic link is not publishable: ' + str(relative))
        if any(part.startswith('.') for part in relative.parts) and relative.as_posix() != '.nojekyll':
            raise ValueError('hidden content is not publishable: ' + str(relative))
        if any(Path(part).stem.lower() in SECRET_NAMES for part in relative.parts):
            raise ValueError('credential-like filename is not publishable: ' + str(relative))
        if stat.S_ISDIR(info.st_mode):
            continue
        if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
            raise ValueError('only ordinary files with no hard links are publishable: ' + str(relative))
        if path.name == 'COMPLETED':
            if info.st_size != len(COMPLETION_MARKER) or path.read_bytes() != COMPLETION_MARKER:
                raise ValueError('COMPLETED must be the exact static-gallery marker: ' + str(relative))
        elif path.suffix.lower() not in EXTENSIONS and relative.as_posix() != '.nojekyll':
            raise ValueError('unsupported static asset; keep raw outputs outside gallery: ' + str(relative))
        if info.st_size > FILE_LIMIT:
            raise ValueError('asset exceeds 90 MB: ' + str(relative))
        total += info.st_size
        if total >= TOTAL_LIMIT:
            raise ValueError('gallery must remain smaller than 1 GB')
        tail = b''
        with path.open('rb') as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b''):
                sample = tail + block
                if SECRET_PATTERN.search(sample):
                    raise ValueError('possible credential content; review locally: ' + str(relative))
                tail = sample[-256:]
        count += 1
    return count, total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('gallery', nargs='?', default='gallery')
    args = parser.parse_args()
    try:
        count, total = check_gallery(args.gallery)
    except (OSError, ValueError) as error:
        parser.exit(1, 'Gallery check failed: ' + str(error) + '\n')
    print('Gallery check passed: {} files, {:.2f} MB. Review content before committing.'.format(count, total / 1e6))


if __name__ == '__main__':
    main()
