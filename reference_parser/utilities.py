
def create_reverse_complement(string):
    new_list = []
    for x in string:
        if x == 'A':
            new_list.append('T')
        if x == 'C':
            new_list.append('G')
        if x == 'G':
            new_list.append('C')
        if x == 'T':
            new_list.append('A')
    new_string = ''.join(new_list)
    return new_string[::-1]

