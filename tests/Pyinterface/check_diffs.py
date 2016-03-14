from filecmp import cmp

if __name__ == "__main__":
    path = './'
    prefix = ['idoff.','idon.']
    suffix = '.dat'
    compare = ['f','py']

    for p in prefix:
        file0, file1 = [path+p+c+suffix for c in compare]
        result = cmp(file0, file1)
        if not result:
            print ': {0}, {1} DO NOT match!!'.format(file0.strip('./'), file1.strip('./'))
        else:
            print ': {0}, {1} are the same'.format(file0.strip('./'), file1.strip('./'))

    path = './'
    prefix = 'fort.'
    suffix = ''
    compare = [('69','71'), ('70','72')]

    for c in compare:
        file0, file1 = [path+prefix+cc for cc in c]
        result = cmp(file0, file1)
        if not result:
            print '| {0}, {1} DO NOT match!!'.format(file0.strip('./'), file1.strip('./'))
        else:
            print '| {0}, {1} are the same'.format(file0.strip('./'), file1.strip('./'))
