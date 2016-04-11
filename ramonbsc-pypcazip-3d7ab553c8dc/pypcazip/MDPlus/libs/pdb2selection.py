def pdb2selection(pdbfile):
    '''
    A little utility function to convert 'mask' pdb files into
    MDAnalysis-style selection strings. Basically it reads the second column
    (atom number) and uses it to construct 'bynum' selections. Runs
    of consecutive numbers are expressed in start:stop form.
    '''

    sel = ''
    i = 0
    j = 0
    with open(pdbfile, 'r') as f:
        for line in f:
            if line.find('ATOM') == 0 or line.find('HETATM') == 0:
                k = int(line.split()[1])
                # the next line catches the initialization process:
                if i == 0:
                    i = k
                    j = k - 1
                # are we in a run of consecutive numbers?:
                if k == j + 1:
                    j = k
                else:
                    # time to write out another selection:
                    sel = sel + ' bynum {0}:{1} or'.format(i, j)
                    i = k
                    j = k
                # end-of-file reached. Make sure last selection is included:
    if i > 0 and j > 0:
        sel = sel + ' bynum {0}:{1} or'.format(i, j)
    if len(sel) > 3:
        # remove the trailing ' or':
        sel = sel[:-3]
    return sel

