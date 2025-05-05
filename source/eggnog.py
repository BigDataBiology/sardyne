def extract_og(ogs):
    has_2 = None
    for c in ogs.split(','):
        og,taxon = c.split('@')
        if taxon == '1':
            return og
        elif taxon == '2':
            has_2 = og
    if has_2:
        return has_2
    if ogs: return ogs.split(',')[0].split('@')[0]
