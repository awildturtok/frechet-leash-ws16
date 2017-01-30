def unique(elements, key):
    test = set()
    out = []

    for element in elements:
        k = key(element)
        if k not in test:
            set.add(k)
            out.append(element)

    return out