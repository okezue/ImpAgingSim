"""Download and verify the project's Zenodo archive (concept DOI 10.5281/zenodo.20499120).

    python scripts/fetch_zenodo.py --dest output/zenodo                 # all files
    python scripts/fetch_zenodo.py --only fixed_density_campaign.tar    # one file
    python scripts/fetch_zenodo.py --extract 'seed_extension/k0.5_'     # extract matching members

The concept DOI always resolves to the newest version; the record id is looked up through
the Zenodo API so the checksums come from the same version that is downloaded.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import tarfile
import urllib.parse
import urllib.request

CONCEPT_DOI = "10.5281/zenodo.20499120"
API = "https://zenodo.org/api/records/"


def latest_record() -> dict:
    query = urllib.parse.quote(f'conceptdoi:"{CONCEPT_DOI}"')
    with urllib.request.urlopen(f"{API}?q={query}") as response:
        hits = json.load(response)["hits"]["hits"]
    if not hits:
        raise RuntimeError(f"no Zenodo record found for concept DOI {CONCEPT_DOI}")
    return hits[0]


def md5_file(path: str) -> str:
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(8 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, dest: str) -> None:
    tmp = dest + ".part"
    with urllib.request.urlopen(url) as response, open(tmp, "wb") as out:
        total = int(response.headers.get("Content-Length") or 0)
        done = 0
        for chunk in iter(lambda: response.read(8 << 20), b""):
            out.write(chunk)
            done += len(chunk)
            if total:
                print(f"\r  {os.path.basename(dest)}: {done/1e6:,.0f}/{total/1e6:,.0f} MB", end="", flush=True)
    print()
    os.replace(tmp, dest)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--dest", default="output/zenodo")
    parser.add_argument("--only", nargs="*", default=None, help="file names to fetch (default: all)")
    parser.add_argument("--extract", default=None, help="substring; extract matching tar members under DEST/extracted")
    args = parser.parse_args(argv)
    record = latest_record()
    print(f"record {record['id']} version {record['metadata'].get('version')} doi {record.get('doi')}")
    os.makedirs(args.dest, exist_ok=True)
    for entry in record["files"]:
        name = entry["key"]
        if args.only is not None and name not in args.only:
            continue
        path = os.path.join(args.dest, name)
        expected = entry["checksum"].split(":", 1)[-1]
        if os.path.exists(path) and md5_file(path) == expected:
            print(f"  {name}: present, checksum ok")
        else:
            download(entry["links"]["self"], path)
            actual = md5_file(path)
            if actual != expected:
                raise RuntimeError(f"checksum mismatch for {name}: {actual} != {expected}")
            print(f"  {name}: downloaded, checksum ok")
        if args.extract and name.endswith(".tar"):
            target = os.path.join(args.dest, "extracted")
            with tarfile.open(path) as archive:
                members = [m for m in archive.getmembers() if args.extract in m.name]
                archive.extractall(target, members=members)
            print(f"  extracted {len(members)} members matching {args.extract!r} into {target}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
