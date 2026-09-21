#!/usr/bin/env python3
"""Conservative local release check; findings report locations, never identifiers.

Run against a checkout with Git LFS files downloaded. This detects known identifier
formats and saved state in cohort-related notebooks, not every possible disclosure.
It does not establish permission to redistribute a dataset or inspect image pixels.
No pickled models are executed. Gzip and ZIP members are scanned after decompression.
"""
import argparse
import csv
import gzip
import io
import json
from pathlib import Path
import re
import sys
import zipfile


ID = re.compile(rb"(?<![A-Za-z0-9])(?:SSC\d{4,}|SP\d{6,}|AU\d{5,}|MSSNG\d{3,}(?:-\d+)?|REACH\d{4,}|SS\d{6,}|SF\d{6,}|SJD[_-]\d+(?:\.\d+)?|\d{1,2}-\d{4}-\d{3}[A-Z]?)(?![A-Za-z0-9])")
CONTEXT = re.compile(r"SSC|SPARK|MSSNG|PsychENCODE|PsychENCONDE|BrainSpan|(?:genome|sample|family|subject|individual|donor)[_ .-]id|proband|pedigree|\.ped\b|VCF_files|patient", re.I)
NUMERIC_ID = re.compile(r"(?:family|sample|genome|proband|donor|subject|individual)[_ .]?(?:id)?\s*(?:==|:=|=|:)\s*[\"']\d[\w.-]*[\"']", re.I)
PERSONAL_PATH = re.compile(r"/Users/|/Volumes/|OneDrive-[^/]+/|Dropbox/")
SECRET = re.compile(rb"(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{40,}|-----BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY-----)")
MOUSE_TABLES = {
    'DataAnalysis/Input/MicEvents.mm10.neural.tab',
    'microexon-code/resources/inputs/wt/MicEvents.mm10.tab.gz',
}


def strings(obj):
    if isinstance(obj, str):
        yield obj
    elif isinstance(obj, list):
        for value in obj:
            yield from strings(value)
    elif isinstance(obj, dict):
        for key, value in obj.items():
            if key.startswith(('image/', 'audio/', 'video/')) or key == 'application/pdf':
                continue
            yield from strings(value)


def scan(root):
    findings = []
    counts = {'files': 0, 'notebooks': 0, 'gzip_files': 0, 'zip_members': 0}

    def flag(path, rule):
        safe_path = ID.sub(b'[REDACTED]', str(path).encode()).decode()
        findings.append({'path': safe_path, 'rule': rule})

    def check_stream(handle, location):
        tail = b''
        id_hit = secret_hit = False
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            data = tail + chunk
            id_hit = id_hit or bool(ID.search(data))
            secret_hit = secret_hit or bool(SECRET.search(data))
            tail = data[-4096:]
        if id_hit:
            flag(location, 'identifier-pattern')
        if secret_hit:
            flag(location, 'credential-pattern')

    for path in sorted(root.rglob('*')):
        rel = path.relative_to(root)
        if '.git' in rel.parts:
            continue
        if path.is_symlink():
            flag(rel, 'symlink-needs-review')
            continue
        if not path.is_file():
            continue
        counts['files'] += 1
        if any(x in rel.parts for x in ('__pycache__', '.ipynb_checkpoints', '.pytest_cache')):
            flag(rel, 'cached-artifact')
        if ID.search(str(rel).encode()):
            flag(rel, 'identifier-in-filename')
        with path.open('rb') as handle:
            if handle.read(150).startswith(b'version https://git-lfs.github.com/spec/v1'):
                flag(rel, 'unhydrated-lfs-pointer')
                continue
        if path.suffix == '.ipynb':
            counts['notebooks'] += 1
            try:
                doc = json.loads(path.read_text())
                text = '\n'.join(strings(doc))
                check_stream(io.BytesIO(text.encode()), rel)
                if PERSONAL_PATH.search(text):
                    flag(rel, 'personal-filesystem-path')
                if NUMERIC_ID.search(text):
                    flag(rel, 'hardcoded-numeric-participant-filter')
                if 'widgets' in doc.get('metadata', {}):
                    flag(rel, 'saved-widget-state')
                if CONTEXT.search(text):
                    for i, cell in enumerate(doc['cells']):
                        if cell.get('outputs') or cell.get('attachments') or cell.get('execution_count') is not None:
                            flag(f'{rel}:cell-{i}', 'saved-cohort-notebook-state')
            except (ValueError, KeyError) as error:
                flag(rel, 'invalid-notebook-' + type(error).__name__)
        elif str(rel) in MOUSE_TABLES:
            # AU-prefixed mouse gene symbols are not human participant identifiers.
            # Exempt only that column in these two known mouse event tables.
            opener = gzip.open if path.suffix == '.gz' else open
            with opener(path, 'rt') as handle:
                for row in csv.DictReader(handle, delimiter='\t'):
                    other_fields = '\t'.join(v or '' for k, v in row.items() if k != 'geneName')
                    check_stream(io.BytesIO(other_fields.encode()), rel)
            counts['gzip_files'] += int(path.suffix == '.gz')
        elif path.suffix == '.gz':
            counts['gzip_files'] += 1
            with gzip.open(path, 'rb') as handle:
                check_stream(handle, rel)
        elif zipfile.is_zipfile(path):
            with zipfile.ZipFile(path) as archive:
                for member in archive.infolist():
                    counts['zip_members'] += 1
                    with archive.open(member) as handle:
                        check_stream(handle, f'{rel}:archive-member')
        else:
            with path.open('rb') as handle:
                check_stream(handle, rel)
    return {'counts': counts, 'findings': findings}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', nargs='?', type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    if not args.directory.is_dir():
        parser.error('directory does not exist')
    report = scan(args.directory.resolve())
    print(json.dumps(report, indent=2))
    return 1 if report['findings'] else 0


if __name__ == '__main__':
    sys.exit(main())
