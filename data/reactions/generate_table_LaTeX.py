f_in = open('ireactions.csv','r')
f_out = open('ireactions_table.tex','w')

t = f_in.readline().split(',')
n_reactions = int(t[0])
f_in.readline()
f_in.readline()

f_out.write('\\documentclass{article}\n')
f_out.write('\\usepackage{fullpage}\n')
f_out.write('\\usepackage{longtable}\n')
f_out.write('\\usepackage[landscape]{geometry}\n')
f_out.write('\\begin{document}\n')
f_out.write('\\begin{longtable}{l l l l l l}\n')
f_out.write('\# & Reaction & $k$ & Reference & Notes\\\\\n')
for i_reactions in range(n_reactions):
    t = f_in.readline().split(',')
    for i_spe in range(1,6):
    # convert to standard symbols
        if t[i_spe] == 'N2D':
            t[i_spe] = 'N$^2$D'
        elif t[i_spe] == 'O1D':
            t[i_spe] = 'O$^1$D'
        elif t[i_spe] == 'E':
            t[i_spe] = 'e$^-$'
        else:
            if t[i_spe].find('P') != -1:
                t[i_spe] = t[i_spe].replace('P','$^+$')
            if t[i_spe].find('2') != -1:
                t[i_spe] = t[i_spe].replace('2','$_2$')
            if t[i_spe].find('3') != -1:
                t[i_spe] = t[i_spe].replace('3','$_3$')
            if t[i_spe].find('$$') != -1:
                t[i_spe] = t[i_spe].replace('$$','')
    f_out.write(str(i_reactions+1) + ' & ')
    # output chemical formula
    if t[2] == '-':
        f_out.write(t[1] + ' $\\rightarrow$ ' + t[3] + ' & ')
    elif t[4] == '-':
        f_out.write(t[1] + ' + ' + t[2] + ' $\\rightarrow$ ' + t[3] + ' & ')
    elif t[5] == '-':
        f_out.write(t[1] + ' + ' + t[2] + ' $\\rightarrow$ ' + t[3] + ' + ' + t[4] + ' & ')
    else:
        f_out.write(t[1] + ' + ' + t[2] + ' $\\rightarrow$ ' + t[3] + ' + ' + t[4] + ' + ' + t[5] + ' & ')
    
    # output reaction type
    #f_out.write(t[6] + ' & ')
    itype = int(t[6])
    
    # output rate equation
    for i_rc in range(7,17):
        if t[i_rc].find('E') != -1:
            t[i_rc] = t[i_rc][0:4]+'\\times 10^{' + str(int(t[i_rc][5:8])) + '}'
    if itype == 1: # Unimolecular
        f_out.write('$' + t[7] + '$ & ')
    elif itype == 2: # Bimolecular
        f_out.write('$' + t[7])
        if float(t[8]) != 0:
            f_out.write('T^{' + t[8] + '}')
        if float(t[9]) != 0:
            f_out.write('e^{' + t[9] + '/T}')
        f_out.write('$ & ')
    elif itype == 3: # Association
        f_out.write('\parbox[t]{7 cm}{')
        if float(t[16]) == 0:
            f_out.write('$\\frac{k_0[M]k_\\infty}{k_0[M]+k_\\infty}')
        else:
            f_out.write('$\\frac{k_0[M]k_\\infty}{k_0[M]+k_\\infty}F')
        f_out.write('$,\\\\')
        f_out.write('$k_\\infty=' + t[7])
        if float(t[8]) != 0:
            f_out.write('T^{' + t[8] + '}')
        if float(t[9]) != 0:
            f_out.write('e^{' + t[9] + '/T}')
        f_out.write('$,\\\\ $k_0=' + t[10])
        if float(t[11]) != 0:
            f_out.write('T^{' + t[11] + '}')
        if float(t[12]) != 0:
            f_out.write('e^{' + t[12] + '/T}')
        f_out.write('$')
        if float(t[16]) != 0:
            f_out.write(',\\\\ $\\log F = \\log F_c\\left(1+(\\frac{\\log(k_0[M]k_\\infty+C)}{(N-0.14(\\log(k_0[M]k_\\infty+C)))})^2\\right)^{-1}$,\\\\')
            f_out.write('$F_c=' + t[16] + ', N = 0.75-1.27\\log F_c, $\\\\ $C = -0.4-0.67 \\log F_c$')
        f_out.write('} &')
    elif itype == 4: # Association# Association & Radiative Association
        f_out.write('\parbox[t]{7 cm}{')
        if float(t[16]) == 0:
            f_out.write('$\\frac{k_0[M]k_\\infty}{k_0[M]+k_\\infty}')
        else:
            f_out.write('$\\frac{k_0[M]k_\\infty}{k_0[M]+k_\\infty}F')
        if t[13] != '0.00\\times 10^{0}':
            f_out.write('+k_R')
        f_out.write('$,\\\\')
        f_out.write('$k_\\infty=' + t[7])
        if float(t[8]) != 0:
            f_out.write('T^{' + t[8] + '}')
        if float(t[9]) != 0:
            f_out.write('e^{' + t[9] + '/T}')
        f_out.write('$,\\\\ $k_0=' + t[10])
        if float(t[11]) != 0:
            f_out.write('T^{' + t[11] + '}')
        if float(t[12]) != 0:
            f_out.write('e^{' + t[12] + '/T}')
        f_out.write('$')
        if t[13] != '0.00\\times 10^{0}':
            f_out.write(',\\\\ $k_R=' + str(t[13]))
            if float(t[14]) != 0:
                f_out.write('T^{' + t[14] + '}')
            if float(t[15]) != 0:
                f_out.write('e^{' + t[15] + '/T}')
            f_out.write('$ ')
        if float(t[16]) != 0:
            f_out.write(',\\\\ $\\log F = \\log F_c\\left(1+(\\frac{\\log(k_0[M]k_\\infty+C)}{(N-0.14(\\log(k_0[M]k_\\infty+C)))})^2\\right)^{-1}$,\\\\')
            f_out.write('$F_c=' + t[16] + ', N = 0.75-1.27\\log F_c, $\\\\ $C = -0.4-0.67 \\log F_c$')
        f_out.write('} &')
    elif itype == 5: # Association (Sander's Formula)
        f_out.write('\parbox[t]{7 cm}{')
        f_out.write('$\\frac{k_0[M]k_\\infty}{k_0[M]+k_\\infty}0.6^\\alpha$, \\\\')
        f_out.write('$k_\\infty = ' + t[7])
        if float(t[8]) != 0:
            f_out.write('T^{' + t[8] + '}')
        if float(t[9]) != 0:
            f_out.write('e^{' + t[9] + '/T}, ')
        f_out.write('$, \\\\$k_0=' + t[10])
        if float(t[11]) != 0:
            f_out.write('T^{' + t[11] + '}')
        if float(t[12]) != 0:
            f_out.write('e^{' + t[12] + '/T')
        f_out.write('$, \\\\ $\\alpha = 1 + \\left(\\log_{10}\\frac{k_0[M]}{k_\\infty}\\right)^2$} &')
    elif itype == -1: # Unimolecular Reaction
        f_out.write('$' + t[7])
        if float(t[8]) != 0:
            f_out.write('T^{' + t[8] + '}')
        if float(t[9]) != 0:
            f_out.write('e^{' + t[9] + '/T}')
        f_out.write('$ & ')
    elif itype == -2: # Two-Body Reaction
        f_out.write('$' + t[7])
        if float(t[8]) != 0:
            f_out.write('T^{' + t[8] + '}')
        if float(t[9]) != 0:
            f_out.write('e^{' + t[9] + '/T}')
        f_out.write('$ & ')
    elif itype == -4: # Electron Recombination
        f_out.write('$' + t[7])
        if float(t[8]) != 0:
            f_out.write('T^{' + t[8] + '}')
        if float(t[9]) != 0:
            f_out.write('e^{' + t[9] + '/T}')
        f_out.write('$ & ')

    # output ancillary info
    f_out.write(t[18])
    #if t[19] != '\n':
    #    f_out.write(' & ' + t[19].replace('\\n',''))
    f_out.write(' \\\\\n')
    

f_out.write('\\end{longtable}')
f_out.write('\\end{document}')
f_in.close()

f_out.close()
