#!/usr/bin/python3
###############################################################################


def remove_label_after_do(str):
    splitted_loc = str.split()
    if splitted_loc[0] == 'do':
        del splitted_loc[1]
    else:
        del splitted_loc[2]
    return ' '.join(splitted_loc)


###############################################################################
def detect_is_fixed_form(file_lines):
    for file_line in file_lines:
        file_line = file_line.lower()
        if file_line[0] == '!' or file_line[0] == 'c':
            continue
        if len(file_line[:-1]) > 72:
            print(file_line)
            return False
        lablel_loc = file_line[0:5].strip()
        if len(lablel_loc) == 0:
            continue
        try:
            int(lablel_loc)
        except:
            print(file_line)
            return False
    return True


###############################################################################
def get_goto_labels(lines):
    goto_label_stack = []
    # direct goto (or go to) statements
    for line in lines:
        line = line.lower()
        if len(line) != 0 and (line[0] == 'c' or line[0] == '!'):
            continue
        if not (('go' in line and 'to' in line) or ('goto' in line)):
            continue
        splited = line.split()
        try:
            label = int(splited[-1])
            goto_label_stack.append(label)
        except:
            continue
    # goto from if
    for line in lines:
        line = line.lower()
        if len(line) != 0 and (line[0] == 'c' or line[0] == '!'):
            continue
        if 'if' not in line:
            continue
        last_bracket = line.rfind(')')+1
        label_string = line[last_bracket:]
        if label_string.count(',') != 2:
            continue
        labels = label_string[:-1].split(',')
        for label in labels:
            goto_label_stack.append(int(label))

    return goto_label_stack


###############################################################################
def reformat_do(lines):
    out_lines = ''
    indent = 0
    do_label_stack = []
    goto_label_stack = get_goto_labels(lines)
    for line in lines:
        isMod = False
        cont_line = ''
        line = line.lower()
        if len(line) != 0 and (line[0] == 'c' or line[0] == '!'):
            line = '!'+line[1:]
            out_lines += line
            continue
        splited = line.split()
        while True:
            if len(splited)==0 : break
            if len(splited)==1 and splited[0] != 'do': break
            if splited[0] != 'do' and splited[1] != 'do': break
            if splited[0] == 'do': label = splited[1]
            else: label = splited[2]
            try:
                label = int(label)
            except:
                break
            do_label_stack.append(label)
            line = lead+'  '*indent+remove_label_after_do(line)+'\n'
            indent += 1
            isMod = True
            break
        label = line[:ii-1]
        try:
            label = int(label)
            if label not in goto_label_stack:
                if line.split()[-1] == "continue":
                    cont_line = line
                    line = ''
                    isMod = True
                elif len(do_label_stack) !=0 and do_label_stack[-1] == label:
                    line = lead+'  '*indent+line[ii:]
                    isMod = True
        except:
            pass
        while len(do_label_stack) !=0 and do_label_stack[-1] == label:
            del do_label_stack[-1]
            indent -= 1
            isMod = True
            line = line+lead+'  '*indent+'end do\n'
            if len(do_label_stack) == 0 or do_label_stack[-1] != label: break
        # if len(cont_line) !=0 and label not in do_label_stack:
        #     line += cont_line
        if not isMod and len(line)>ii-1: line = line[:ii]+'  '*indent+line[ii:]
        if len(line[:-1]) > 72:
            line = line[:72]+'\n     &'+line[72:]
        # print(line[:-1])
        out_lines += line
    return out_lines


###############################################################################
###############################################################################
fname = "from.txt"

isFixedForm = True
with open(fname, 'r') as fin:
    lines = fin.readlines()
    isFixedForm = detect_is_fixed_form(lines)
    print("Is fixed form:", isFixedForm)



ext = '.f90'
lead = ''
ii = 0
if isFixedForm:
    ext = '.f'
    lead = '      '
    ii = 6

with open(fname,'r') as fin:
    with open (fname[:-4]+"_1"+ext,'w') as fout:
        lines = fin.readlines()
        out_lines = reformat_do(lines)
        fout.write(out_lines)
