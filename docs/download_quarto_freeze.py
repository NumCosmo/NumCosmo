#!/usr/bin/env python3

import json
import os
import sys
import zipfile
import requests

REPO = "NumCosmo/NumCosmo"
SHA = os.environ["READTHEDOCS_GIT_COMMIT_HASH"]

artifact_name = f"quarto-freeze-{SHA}"

print(f"Looking for {artifact_name}")

url = f"https://api.github.com/repos/{REPO}/actions/artifacts"

r = requests.get(url)
r.raise_for_status()

artifacts = r.json()["artifacts"]

artifact = None
for a in artifacts:
    if a["name"] == artifact_name and not a["expired"]:
        artifact = a
        break

if artifact is None:
    print(f"Artifact not found: {artifact_name}")
    sys.exit(1)

download_url = artifact["archive_download_url"]

print(f"Downloading {artifact_name}")

r = requests.get(
    download_url,
    headers={"Accept": "application/vnd.github+json"},
)
r.raise_for_status()

with open("freeze.zip", "wb") as f:
    f.write(r.content)

target = "build/docs/.quarto/_freeze"

os.makedirs(target, exist_ok=True)

with zipfile.ZipFile("freeze.zip") as z:
    z.extractall(target)

print(f"Extracted into {target}")
