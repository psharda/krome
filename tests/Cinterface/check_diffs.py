from filecmp import cmp

if __name__ == "__main__":

    path = './'
    prefix = 'fort.'
    suffix = ''
    compare = [('69','71'), ('70','72'), ('66','77'), ('67','78')]

    for c in compare:
        file0, file1 = [path+prefix+cc for cc in c]
        result = cmp(file0, file1)
        if not result:
            print '| {0}, {1} DO NOT match!!'.format(file0.strip('./'), file1.strip('./'))
        else:
            print '| {0}, {1} are the same'.format(file0.strip('./'), file1.strip('./'))
