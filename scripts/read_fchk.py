import re
def fchk_fields(path):
    lines = path.read_text().splitlines()
    fields = {}
    i = 0
    while i < len(lines):
        match = re.match(r'^(.+?)\s{2,}([IR])\s+(.*)$', lines[i])
        i += 1
        if not match: continue
        key, kind, payload = match.groups()
        convert = int if kind == 'I' else lambda value: float(value.replace('D','E'))
        if payload.startswith('N='):
            n = int(payload[2:])
            values = []
            while len(values) < n:
                values.extend(convert(x) for x in lines[i].split())
                i += 1
            if len(values) != n: raise ValueError(f'Malformed {key}')
            fields[key] = values
        else:
            fields[key] = convert(payload.strip())
    return fields
