def linePrepend(filename,line):
    '''
    simple script to prepend a line at the being of all files in a folder
    usually will be a line containing column headers
    '''
    with open(filename, 'r+') as file:
        content = file.read()
        file.seek(0,0)
        file.write(line + '\n' + content)
        file.close()